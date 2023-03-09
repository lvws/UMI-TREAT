#!/usr/bin/env python
import os,sys,subprocess
import argparse
from subprocess import PIPE

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input",help="name sorted bam",required=True)
    parser.add_argument("-p","--prefix",help="out file prefix")
    parser.add_argument("-t","--threads",help="threads run",default=8)
    parser.add_argument("--pigz",help="pizg path",)
    parser.add_argument("-o","--outfold",help="output file fold.",required=True)
    args = parser.parse_args()
    if not args.prefix:
        args.prefix = os.path.basename(args.input).split('.')[0]
    if not args.pigz:
        script_path = os.path.dirname(os.path.abspath(sys.argv[0]))
        args.pigz = "{script_path}/tools/pigz".format(**locals())
    return args

def runUMItreat(args):
    script_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    bam = args.input
    prefix = args.prefix
    threads = args.threads
    outfold = args.outfold
    cmd = "{script_path}/target/release/bamtest {bam} -p {prefix} -t {threads} -o {outfold}".format(**locals())
    run = subprocess.Popen(cmd,stdout=PIPE,stderr=PIPE,shell=True)
    for i in run.stdout:
        print(i)
    for i in run.stderr:
        print(i)
    
def combine(args):
    pigz = args.pigz
    prefix = args.prefix
    threads = args.threads
    outfold = args.outfold
    scs_r1_lst = []
    scs_r2_lst = []
    dcs_r1_lst = []
    dcs_r2_lst = []
    for i in range(1,int(args.threads)+1):
        scs_r1_lst.append("{0}/{1}_r1-scs-{2}.fastq".format(args.outfold,args.prefix,i))
        scs_r2_lst.append("{0}/{1}_r2-scs-{2}.fastq".format(args.outfold,args.prefix,i))
        dcs_r1_lst.append("{0}/{1}_r1-dcs-{2}.fastq".format(args.outfold,args.prefix,i))
        dcs_r2_lst.append("{0}/{1}_r2-dcs-{2}.fastq".format(args.outfold,args.prefix,i))

    scs_r1_files = " ".join(scs_r1_lst)
    scs_r2_files = " ".join(scs_r2_lst)
    dcs_r1_files = " ".join(dcs_r1_lst)
    dcs_r2_files = " ".join(dcs_r2_lst)

    cmd1 = "cat {scs_r1_files} | {pigz} > {outfold}/{prefix}_scs_R1.fastq.gz".format(**locals())
    cmd2 = "cat {scs_r2_files} | {pigz} > {outfold}/{prefix}_scs_R2.fastq.gz".format(**locals())
    cmd3 = "cat {dcs_r1_files} | {pigz} > {outfold}/{prefix}_dcs_R1.fastq.gz".format(**locals())
    cmd4 = "cat {dcs_r2_files} | {pigz} > {outfold}/{prefix}_dcs_R2.fastq.gz".format(**locals())

    p1 = subprocess.Popen(cmd1,shell=True)
    p2 = subprocess.Popen(cmd2,shell=True)
    p3 = subprocess.Popen(cmd3,shell=True)
    p4 = subprocess.Popen(cmd4,shell=True)
    p1.wait()
    p2.wait()
    p3.wait()
    p4.wait()

    os.system("rm {scs_r1_files} {scs_r2_files} {dcs_r1_files} {dcs_r2_files}".format(**locals()))

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.argv.append('-h')
    #script_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    args = get_args()
    print("START TREAT...")
    runUMItreat(args)
    print("FINISHED TREAT.")
    print("START COMBINE...")
    combine(args)
    print("FINISHED!")