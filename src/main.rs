use rust_htslib::bam::{self,Read, Record};
use std::io::Write;
use std::{error::Error, cmp::min,cmp::max, default};
use std::sync::{Arc,Mutex};
use rustc_hash::FxHashMap;
use ahash::AHashMap;
use std::fs::File;
use std::thread;
use clap::{Command,Arg,ArgAction};

type MyResult<T> = Result<T,Box<dyn Error>>;
type vreads = Vec<ReadI>;
type ureadsDic =  AHashMap<String,vreads>;
type tDic = AHashMap<String,ureadsDic>;
type tVec = Vec<(String,ureadsDic)>;


struct Config{
    input:String,
    outfold:String,
    prefix:String,
    threads:usize,
}

fn getConfig() -> Config{
    let mut matches = Command::new("umiTreat")
        .about("UMI MisMatch 1 ,make scs and dcs")
        .author("lvws; 1129309258@qq.com")
        .version("0.01")
        .arg(Arg::new("INPUT")
            .help("name sorted bam , input to make scs and dcs.")
            .required(true)
        )
        .arg(Arg::new("OUTFOLD")
            .long("outfold")
            .short('o')
            .help("out file store fold")
            .required(true)
        )
        .arg(Arg::new("PREFIX")
            .long("prefix")
            .short('p')
            .help("out file prefix")
            .required(true)
        )
        .arg(Arg::new("THREADS")
            .long("threads")
            .short('t')
            .help("threads to run,default is 8")
            .default_value("8")
            .value_parser(clap::value_parser!(usize))
        )
        .get_matches();

    Config { input: matches.remove_one::<String>("INPUT").unwrap(), 
        outfold: matches.remove_one::<String>("OUTFOLD").unwrap(), 
        prefix: matches.remove_one::<String>("PREFIX").unwrap(), 
        threads: *matches.get_one("THREADS").unwrap(), 
        }
}

fn main() -> MyResult<()>{
    println!("START UMI TREAT!");
    let config = getConfig();
    let threads = config.threads;
    let outfold = config.outfold;
    let prefix = config.prefix;
    let bam = config.input;
    let mut pos_dic:tDic = AHashMap::new();
    read_bam(&bam,&mut pos_dic)?;
    println!("Finish pos_dic!");
    let mut pos_vec = Vec::from_iter(pos_dic);
    println!("Fish trans_dic2vec");
    let num = pos_vec.len()/threads;
    let mut handles = vec![];
    for index in 1..=threads-1{
        //let pos_vec_arc = Arc::new(Mutex::new(pos_vec.split_off(num*(threads-index))));
        let pos_vec_arc = pos_vec.split_off(num*(threads-index));
        println!("Start: thread {}",index);
        let outfile = format!("{outfold}/{prefix}");
        let tmp_thread = thread::spawn(move || {
            parseReadI(pos_vec_arc,outfile, index).unwrap();
        });
        handles.push(tmp_thread);
    }
    //let pos_vec_arc = Arc::new(Mutex::new(pos_vec));
    let outfile = format!("{outfold}/{prefix}");
    let temp_thread = thread::spawn(move||{
        parseReadI(pos_vec, outfile,threads).unwrap();
    });

    handles.push(temp_thread);
    for h in handles{
        h.join().unwrap();
    }
    println!("Finshed!");
    Ok(())
}

/*
记录bam read 信息；
 */
#[derive(Debug)]
pub struct ReadI{
    umi_r1:String,
    umi_r2:String,
    flag_r1:u16,
    flag_r2:u16,
    qname:String,
    seq_r1:Vec<u8>,
    qual_r1:Vec<u8>,
    seq_r2:Vec<u8>,
    qual_r2:Vec<u8>,
    pos_key:String,
    utag:String,
    cons:u8,//记录sscs 的合并数
    dcs:u8,//记录dcs的另一个ssc 数
    //utag:String, //记录UMI的格式，5-6、7-8之类
}


/*
从read1,read2 中获取序列UMI信息，构建ReadI 结构体
 */
pub fn get_read_info1(urad:&Record,urad2:&Record) -> ReadI{    
    let qname = String::from_utf8_lossy(urad.qname()).to_string();
    let flag = urad.flags();
    // 判断逆向
    let mut seq;
    let mut quals;   
    if flag & 16 != 0 {
        seq = reverse_complement(&urad.seq().as_bytes());
        quals = urad.qual().to_vec();
        quals.reverse();
    } else {
        seq = urad.seq().as_bytes();
        quals = urad.qual().to_vec();
    }

    let tid = urad.tid();
    let pos = urad.pos();
    let mtid = urad.mtid();
    let mpos = urad.mpos();
    let insertSize = urad.insert_size();

    let qid:Vec<&str> = qname.split('|').collect();
    let umis:Vec<&str>= qid[0].split('-').collect();
    let mut umi_r1 = "";
    let mut umi_r2 = "";
    if qid[1] == "ab" {
        umi_r1 = umis[0];
        umi_r2 = umis[1];
    } else {
        umi_r1 = umis[1];
        umi_r2 = umis[0];
    }
    let utag = format!("{}-{}",umi_r1.len(),umi_r2.len());
    umi_r1 = &umi_r1[..3];
    umi_r2 = &umi_r2[..3];
    
    
    /*
    计算比对位置信息，作为key， 汇总相同位置的reads
        */
    let mut pos_key = "".to_string();
    if tid == mtid {
        let start = min(pos, mpos);
        let end = start + insertSize.abs();
        pos_key = format!("{tid}:{start}-{end}");
    } else {
        pos_key = format!("{tid}:{pos}-{mtid}:{mpos}");
    }

    // urad2
    //let qname2 = String::from_utf8_lossy(urad2.qname()).to_string();
    let mut seq2;
    let mut qual2 ;
    let flag2 = urad2.flags();

    if flag2 & 16 != 0{
        seq2 = reverse_complement(&urad2.seq().as_bytes());
        qual2 = urad2.qual().to_vec();
        qual2.reverse();
    } else {
        seq2 = urad2.seq().as_bytes();
        qual2 = urad2.qual().to_vec();
    }

    let mut seq_r1 = vec![];
    let mut seq_r2 = vec![];
    let mut qual_r1:Vec<u8> = Vec::new();
    let mut qual_r2:Vec<u8> = Vec::new() ;
    let mut flag_r1= 0;
    let mut flag_r2 = 0;
    if flag & 64 != 0 {
        seq_r1 = seq;
        seq_r2 = seq2;
        qual_r1 = quals;
        qual_r2 = qual2;
        flag_r1 = flag;
        flag_r2 = flag2;
    } else {
        seq_r1 = seq2;
        seq_r2 = seq;
        qual_r1 = qual2;
        qual_r2 = quals;
        flag_r1 = flag2;
        flag_r2 = flag;
    }

    ReadI { umi_r1:umi_r1[..3].to_string(), umi_r2:umi_r2[..3].to_string(),qname, flag_r1,flag_r2,seq_r1, qual_r1, seq_r2, qual_r2,pos_key,utag,cons:1,dcs:0 }

}

/*
读取BAM 文件获取reads 信息。
 */
pub fn read_bam(infile:&str,pos_dic:&mut tDic) -> MyResult<()>{
    let mut bam = bam::Reader::from_path(infile)?;
    let recodes = bam.records(); 
    let mut pair_read = Vec::new();
    for read in recodes{
        let urad1 = read?;
        if urad1.flags() & 256 == 0{
            pair_read.push(urad1);
        }
        if pair_read.len() == 2{
            let readI = get_read_info1(&pair_read[0], &pair_read[1]);
            if let Some(x) = pos_dic.get_mut(readI.pos_key.as_str()) {
                if let Some(y) = x.get_mut(readI.utag.as_str()){
                    y.push(readI);
                } else {
                    let key = readI.utag.clone();
                    let tmp :vreads = vec![readI];
                    x.insert(key, tmp);
                }
            } else {
                let mut tmp :ureadsDic = AHashMap::new();
                let key = readI.pos_key.clone();
                tmp.insert(readI.utag.clone(), vec![readI]);
                pos_dic.insert(key,tmp);
            }
            pair_read.clear();
        }

    }
    Ok(())
}


/*
解析ReadI 构建SCS 和 DCS 
 */
pub fn parseReadI(pos_vec:tVec,outfile:String,index:usize) -> MyResult<()>{
   
    let mut data = pos_vec;
    let vcc = vec!["5-5".to_string(),"6-6".to_string(),"7-7".to_string(),"8-8".to_string()];
    let mut scs_r1 = File::create(format!("{outfile}_r1-scs-{index}.fastq")).unwrap();
    let mut scs_r2 = File::create(format!("{outfile}_r2-scs-{index}.fastq")).unwrap();
    let mut dcs_r1 = File::create(format!("{outfile}_r1-dcs-{index}.fastq")).unwrap();
    let mut dcs_r2 = File::create(format!("{outfile}_r2-dcs-{index}.fastq")).unwrap();
    for (_,udic) in data.iter_mut(){
        // 构建SSCS；同一个比对位置的所有UMI类型都进行了构建；
        // 返回的{umi:[[],[],[]]} 每个子vector 的第一个元素是构建好的sscs的index位置。
        let mut ssc_dic: FxHashMap<String,Vec<usize>> = FxHashMap::default();
        for (umi,vreads) in udic.iter_mut(){
            let scs_index = findScs(&vreads);
            let mut newv = vec![]; // 记录了SSCS,的index
            for i in scs_index{
                if i.len() > 1{                 
                    for j in 1..i.len(){
                        makeScs(vreads, i[0], i[j]);
                    }                    
                }
                newv.push(i[0]);
            }
            ssc_dic.insert(umi.clone(),newv);
        }

        for utag in vcc.iter(){
            if let Some(v) = udic.get_mut(utag){
                //[[1],[2,3],[4,6],5] 第一项是我们想要的reads,长度为1 的没有DCS， 长度为2的，构建DCS。
                let dcs_index = findDcs_cc(&v, &ssc_dic[utag]);
                for cc in dcs_index{
                    if cc.len() == 1{
                        writer(&v[cc[0]], &mut scs_r1, &mut scs_r2)?;
                    } else{
                        makeDcs_cc(v,cc[0],cc[1]);
                        writer(&v[cc[0]], &mut dcs_r1, &mut dcs_r2)?;
                    }
                }
            }
        }

        // UMI-ab
        // 5-6, 6-5
        let umis = vec![["5-6".to_string(),"6-5".to_string()],["5-7".to_string(),"7-5".to_string()],["5-8".to_string(),"8-5".to_string()],
            ["6-7".to_string(),"7-6".to_string()],["7-8".to_string(),"8-7".to_string()]];

        //let atmp = RefCell::new(udic);
        for umi in umis{
            let mut empty1_readI:&Vec<ReadI> = &vec![];
            let mut empty1_index:&Vec<usize> = &vec![];
            let mut empty2_readI:&Vec<ReadI> = &vec![];
            let mut empty2_index:&Vec<usize> = &vec![];

            //let mut binding = atmp.borrow_mut();
            
            if let Some(u1) = udic.get(&umi[0]){
                empty1_readI = u1;
                empty1_index = ssc_dic.get(&umi[0]).unwrap();
            }
            //let mut binding2 = atmp.borrow_mut();
            if let Some(u2) = udic.get(&umi[1]){
                empty2_readI = u2;
                empty2_index = ssc_dic.get(&umi[1]).unwrap();
            }

            let result = findDcs(&empty1_readI, empty1_index, &empty2_readI, empty2_index);
            
            for i in result{
                if i.len() == 1{
                    writer(&empty1_readI[i[0]], &mut scs_r1, & mut scs_r2)?;
                } else if i.len() == 3 {
                    writer(&empty2_readI[i[0]], &mut scs_r1, & mut scs_r2)?;
                } else if i.len() == 2{
                    let rad = makeDcs(&empty1_readI[i[0]],&empty2_readI[i[1]]);
                    writer(&rad, &mut dcs_r1, & mut dcs_r2)?;
                } 
            }
        }
    }
    println!("Thread {index} Start Clean data!");
    //data.clear();
    unsafe{
        data.set_len(0);
    }
    println!("Thread {index} Finished Clean data!");
    Ok(())
}

/*
最终文件读写处理
 */
pub fn writer(rad:&ReadI,r1f:&mut File,r2f:&mut File) -> MyResult<()>{
    r1f.write(format!("{}|{}|{} 1:N:0\n",rad.qname,rad.cons,rad.dcs).to_string().as_bytes())?;
    r2f.write(format!("{}|{}|{} 2:N:0\n",rad.qname,rad.cons,rad.dcs).to_string().as_bytes())?;
    write_seq(rad.flag_r1, &rad.seq_r1, &rad.qual_r1, r1f)?;
    write_seq(rad.flag_r2, &rad.seq_r2, &rad.qual_r2, r2f)?;
    Ok(())
}

pub fn write_seq(flag:u16,seq:&Vec<u8>,qual:&Vec<u8>,file:&mut File) -> MyResult<()>{
    if flag & 16 == 0{
        file.write(seq)?;
        file.write(&[10])?;
        file.write(&qual.iter().map(|x| x + 33).collect::<Vec<u8>>())?;
        file.write(&[10,43,10])?; // "\n+\n"
    } else{      
        file.write(&reverse_complement(&seq))?;
        file.write(&[10])?;
        let mut tmp:Vec<u8> = vec![];
        for i in (0..qual.len()).rev(){
            tmp.push(qual[i]+33);
        }
        file.write(&tmp)?;
        file.write(&[10,43,10])?;
    }
    Ok(())
}

/*
不同长UMI的DCS 处理
 */
pub fn dcs_ab(umi1:&str,umi2:&str,udic:&ureadsDic,ssc_dic:&FxHashMap<String,Vec<usize>>) -> Vec<Vec<usize>>{
    let mut result:Vec<Vec<usize>> = vec![];
    if let Some(u1) = ssc_dic.get(umi1){
        if let Some(u2) = ssc_dic.get(umi2){
            return findDcs(&udic.get(umi1).unwrap(), u1, &udic.get(umi2).unwrap(), u2);
        } else {
            for i in u1{
                result.push(vec![*i]);
            }
            return result;
        }
    } else if let Some(u2) = ssc_dic.get(umi2){
        for i in u2{
            result.push(vec![*i,*i,*i]);
        }
        return result;
    }
    result
}


/*
misMatch 检查
 */
pub fn misMatch(s1:&str,s2:&str) -> u8{
    let mut mis:u8 = 0;
    for i in [..s1.len()]{
        if s1[i] != s2[i]{
            mis += 1;
        }
    }
    mis
}

/*
找到mis < 2 的vector Index , sscs;
输入的是同一位点，相同UMI tag(5-6、5-7...)的ReadI vector;
输出：同一组的index vector. vector 长度为1 的不能构建sscs；
输出： [[0],[1,2],[3]]
 */
pub fn findScs(vs:&Vec<ReadI>) -> Vec<Vec<usize>>{
    let mut dic:FxHashMap<usize,u8> = FxHashMap::default();
    let mut result:Vec<Vec<usize>> = vec![];

    if vs.len() == 0{
        return result;
    } else if vs.len() == 1{
        result.push(vec![0]);
        return result;
    }

    for i in 0..vs.len()-1{
        let mut tmp_v:Vec<usize> = vec![];
        if dic.contains_key(&i){
            continue;
        }
        tmp_v.push(i);
        for j in i+1..vs.len(){
            /* 防止 r1 r2 umi一致的 reads 错误构建SSCS（可能是DCS）
            */               
            if (vs[i].flag_r1 != vs[j].flag_r1) | (vs[i].flag_r2 != vs[j].flag_r2){
                continue;
            }
           
            if misMatch(&vs[i].umi_r1, &vs[j].umi_r1) + misMatch(&vs[i].umi_r2, &vs[j].umi_r2) < 2{
                tmp_v.push(j);
                dic.insert(j, 1);
                
            }
        }
        result.push(tmp_v);
    }
    // 如果最后一个readI 没有在前面加入到别的组，最后把他加进来
    if !dic.contains_key(&(vs.len()-1)){
        result.push(vec![vs.len()-1]);
    }
    result
}

/*
找到mis < 2 的vector Index , dcs;
当UMI r1=r2 时， flag 应该相反；
输出：同一组的index vector. vector 长度为1 的不能构建dcs；
输出： [[0],[1,2],[3]]
 */
pub fn findDcs_cc(vs:&Vec<ReadI>,index:&Vec<usize>) -> Vec<Vec<usize>>{
    let mut dic:FxHashMap<usize,u8> = FxHashMap::default();
    let mut result:Vec<Vec<usize>> = vec![];
    if index.len() == 0{
        return result;
    } else if index.len() == 1{
        result.push(vec![index[0]]);
        return result;
    }
    for i in 0..index.len()-1{
        let mut tmp_v:Vec<usize> = vec![];
        if dic.contains_key(&i){
            continue;
        }
        tmp_v.push(index[i]);
        for j in i+1..index.len(){
            // if(vs[index[i]].flag_r1 != vs[index[j]].flag_r2) | (vs[index[i]].flag_r2 != vs[index[j]].flag_r1){
            //     continue;
            // }
            if misMatch(&vs[index[i]].umi_r1, &vs[index[j]].umi_r2) + misMatch(&vs[index[i]].umi_r2, &vs[index[j]].umi_r1) < 2{
                tmp_v.push(index[j]);
                dic.insert(j, 1);
                break;
            }
        }
        result.push(tmp_v);
    }
    // 如果最后一个readI 没有在前面加入到别的组，最后把他加进来
    if !dic.contains_key(&(index.len()-1)){
        result.push(vec![index[index.len()-1]]);
    }
    result
}

/*
输入： UMI 5-6 的Vec<ReadI>, 5-6 中合并的SCS Index， UMI 6-5 的Vec<ReadI>, 6-5中合并的SCS Index;
输出： [[0],[1,2],[3,3,3]] 
长度1 为 vec1 的元素scs， 长度2 是vec1 的元素 DSC 用第一个元素 ,长度3 是 vec2的元素，scs；
 */
pub fn findDcs(vs1:&Vec<ReadI>,index1:&Vec<usize>,vs2:&Vec<ReadI>,index2:&Vec<usize>) -> Vec<Vec<usize>>{
    let mut dic:FxHashMap<usize,u8> = FxHashMap::default();
    let mut result:Vec<Vec<usize>> = vec![];
    if vs1.len() == 0{
        for i in index2{
            result.push(vec![*i,*i,*i]);
        }
        return result;
    }
    if vs2.len() == 0{
        for i in index1{
            result.push(vec![*i]);
        }
        return result;
    }

    for i in index1{
        let mut tmp = vec![];
        tmp.push(*i);
        for j in index2{
            // if (vs1[*i].flag_r1 != vs2[*j].flag_r2) | (vs1[*i].flag_r2 != vs2[*j].flag_r2){
            //     continue;
            // }
            if misMatch(&vs1[*i].umi_r1, &vs2[*j].umi_r2) + misMatch(&vs1[*i].umi_r2, &vs2[*j].umi_r1) < 2{
                tmp.push(*j);
                dic.insert(*j, 1);
                break;
            }
        }
        result.push(tmp);
    }
    
    for i in index2{
        if !dic.contains_key(i){
            result.push(vec![*i,*i,*i]);
        }
    }

    result

}
/*
sscs 构建，出现不匹配的碱基赋值N，质量将到最低0。
N 对应u8： 78
 */
pub fn makeScs(vrad:&mut Vec<ReadI>,i1:usize,i2:usize){
    //println!("SCS: {}\t{}",vrad[i1].qname,vrad[i2].qname);
    vrad[i1].cons += 1;
    let r1_len = min(vrad[i1].qual_r1.len(),vrad[i2].qual_r1.len());
    let r2_len = min(vrad[i1].qual_r2.len(),vrad[i2].qual_r2.len());
    for i in 0..r1_len{
        if vrad[i1].seq_r1[i] != vrad[i2].seq_r1[i]{
            vrad[i1].seq_r1[i] = 78;
            vrad[i1].qual_r1[i] = 0;
            vrad[i2].seq_r1[i] = 78;
            vrad[i2].qual_r1[i] = 0;
        }
    }
    for i in 0..r2_len{
        if vrad[i1].seq_r2[i] != vrad[i2].seq_r2[i]{
            vrad[i1].seq_r2[i] = 78;
            vrad[i1].qual_r2[i] = 0;
            vrad[i2].seq_r2[i] = 78;
            vrad[i2].qual_r2[i] = 0;
        }
    }

    if vrad[i1].seq_r1.len() < vrad[i2].seq_r1.len() {
        vrad[i1].seq_r1 = vrad[i2].seq_r1.clone();
        vrad[i1].qual_r1 = vrad[i2].qual_r1.clone();
    }
    if vrad[i1].seq_r2.len() < vrad[i2].seq_r2.len() {
        vrad[i1].seq_r2 = vrad[i2].seq_r2.clone();
        vrad[i1].qual_r2 = vrad[i2].qual_r2.clone();
    }

    //ReadI { seq_r1,seq_r2,qual_r1,qual_r2,cons,..*rad1}
}

/*
DCS 构建。
DCS 符合下面的描述；
1. 模板正链r1 的序列应和 模板 r2 的序列互补，BAM中，表现为 flag 相反，reads 相同； 同理r2。
 */
pub fn makeDcs_cc(vrad:&mut Vec<ReadI>,i1:usize,i2:usize){
    vrad[i1].dcs = vrad[i2].cons;
    //println!("DCS: {}\t{}",vrad[i1].qname,vrad[i2].qname);
    let r1_len = min(vrad[i1].seq_r1.len(),vrad[i2].seq_r2.len());
    let r2_len = min(vrad[i1].seq_r2.len(),vrad[i2].seq_r1.len());

    for i in 0..r1_len{
        
        if vrad[i1].seq_r1[i] != vrad[i2].seq_r2[i]{
            vrad[i1].seq_r1[i] = 78;
            vrad[i1].qual_r1[i] = 0;
            vrad[i2].seq_r2[i] = 78;
            vrad[i2].qual_r2[i] = 0;
        }
    }

    for i in 0..r2_len{
        if vrad[i1].seq_r2[i] != vrad[i2].seq_r1[i]{
            vrad[i1].seq_r2[i] = 78;
            vrad[i1].qual_r2[i] = 0;
            vrad[i2].seq_r1[i] = 78;
            vrad[i2].qual_r1[i] = 0;
        }
    }

    if vrad[i1].seq_r1.len() < vrad[i2].seq_r2.len(){
        vrad[i1].seq_r1 = vrad[i2].seq_r2.clone();
        vrad[i1].qual_r1 = vrad[i2].qual_r2.clone();
    }

    if vrad[i1].seq_r2.len() < vrad[i2].seq_r1.len(){
        vrad[i1].seq_r2 = vrad[i2].seq_r1.clone();
        vrad[i1].qual_r2 = vrad[i2].qual_r1.clone();
    }
    
}

//pub fn makeDcs(rad1:&ReadI,rad2:& ReadI){
pub fn makeDcs(rad1:&ReadI,rad2:&ReadI) -> ReadI{
    //println!("DCS: {}\t{}",rad1.qname,rad2.qname);
    let cons = rad1.cons;
    let dcs = rad2.cons;
    let mut seq_r1:Vec<u8>;
    let mut qual_r1:Vec<u8>;
    let mut seq_r2:Vec<u8>;
    let mut qual_r2:Vec<u8>;
    let seq_r1_ref:&Vec<u8>;
    let seq_r2_ref:&Vec<u8>;

    if rad1.seq_r1.len() < rad2.seq_r2.len(){
        seq_r1 = rad2.seq_r2.clone();
        qual_r1 = rad2.qual_r2.clone();
        seq_r1_ref = &rad1.seq_r1;
    } else {
        seq_r1 = rad1.seq_r1.clone();
        qual_r1 = rad1.qual_r1.clone();
        seq_r1_ref = &rad2.seq_r2;
    }

    if rad1.seq_r2.len() < rad2.seq_r1.len(){
        seq_r2 = rad2.seq_r1.clone();
        qual_r2 = rad2.qual_r1.clone();
        seq_r2_ref = &rad1.seq_r2;
    } else {
        seq_r2 = rad1.seq_r2.clone();
        qual_r2 = rad1.qual_r2.clone();
        seq_r2_ref = &rad2.seq_r1;
    }

    for i in 0..seq_r1_ref.len(){
        if seq_r1[i] != seq_r1_ref[i]{
            seq_r1[i] = 78;
            qual_r1[i] = 0;
        }
    }

    for i in 0..seq_r2_ref.len(){
        if seq_r2[i] != seq_r2_ref[i]{
            seq_r2[i] = 78;
            qual_r2[i] = 0;
        }
    }

    ReadI { cons,dcs,seq_r1,seq_r2,qual_r1,qual_r2,umi_r1:rad1.umi_r1.clone(),umi_r2:rad1.umi_r2.clone(),
        flag_r1:rad1.flag_r1,flag_r2:rad1.flag_r2,qname:rad1.qname.clone(),pos_key:rad1.pos_key.clone(),
        utag:rad1.utag.clone() }
    

}

/*
反向互补：
T: 84
A: 65
C: 67
G: 71
N: 78
 */
pub fn reverse_complement(seq:&Vec<u8>) -> Vec<u8>{
    let mut result= vec![];
    for i in (0..seq.len()).rev(){
        match seq[i] {
            84 => result.push(65),
            65 => result.push(84),
            67 => result.push(71),
            71 => result.push(67),
            _ => result.push(78),
        }
    }
    result
}