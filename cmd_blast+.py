#!/usr/bin/python3

# cmd_blast+.py

# Hsin-Ying Chang <hyhazelchang@gmail.com>
# v1 2023/06/10

# Usage: python3 /home/xinchang/pyscript_xin/cmd_blast+.py --blast_dir=/usr/local/ncbi-blast-2.11.0+/bin/blastp --in_dir=/scratch/xinchang/cyano11/cyano11.21/query/ --out_dir=/scratch/xinchang/cyano11/cyano11.21/blast/ --sh_dir=/scratch/xinchang/cyano11/cyano11.21/sh/blast/ --in_file_ext=fasta --out_file_ext=xml --opt="-db /scratch/xinchang/cyano11/cyano11.21/db/faa/cyano11.21 -task blastp -evalue 1e-15 -outfmt 5" --n_job=20


import argparse
import os
import glob

def main():
    parser = argparse.ArgumentParser(
        description=("Make the shscript for blast execution."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--blast_dir",
                        default=None,
                        type=str,
                        help="BLAST file directory. Please provide absolute path.")
    parser.add_argument("--in_dir",
                        default=None,
                        type=str,
                        help="Directory containing input files (the file should have fasta extension). Please provide absolute path.")
    parser.add_argument("--out_dir",
                        default=None,
                        type=str,
                        help="Output directory. Please provide absolute path.")
    parser.add_argument("--sh_dir",
                        default=None,
                        type=str,
                        help="Directory for shell scripts. Please provide absolute path.")
    parser.add_argument("--in_file_ext",
                        default="fasta",
                        type=str)
    parser.add_argument("--out_file_ext",
                        default="xml",
                        type=str)
    parser.add_argument("--opt",
                        default=None,
                        type=str)    
    parser.add_argument("--n_job",
                        default=1,
                        type=int)

    args = parser.parse_args()
    blast_dir = args.blast_dir
    in_dir = args.in_dir
    out_dir = args.out_dir
    sh_dir = args.sh_dir
    in_file_ext = args.in_file_ext
    out_file_ext = args.out_file_ext
    opt = args.opt
    n_job = args.n_job

    # Find the sequence files
    seqfiles = glob.glob(in_dir + "*." + in_file_ext)

    # Make output directory
    if not os.path.exists(out_dir):
        os.system("mkdir -p " + out_dir)

    # Get file names
    count = 0
    blast_cmd = []
    for seq in seqfiles:
        count += 1
        seq_name = os.path.basename(seq).split(".")[0] + "_" + str(count)
        blast_cmd.append(blast_dir + " -query " + seq + " -out " + out_dir + seq_name + "." + out_file_ext + " " + opt)

    # print out job scripts
    os.system("mkdir -p " + sh_dir)
    quo = int(count / n_job)
    mod = int(count % n_job)
    cmd_num = 0
    for n in range(n_job):
        job = open(sh_dir + "job" + str(n+1) + ".sh" , "w")
        if n + 1 <= mod:
            for num in range(quo + 1):
                job.write(blast_cmd[cmd_num] + "\n")
                cmd_num += 1
            job.close()
        else:
            for num in range(quo):
                job.write(blast_cmd[cmd_num] + "\n")
                cmd_num += 1
            job.close()

if __name__ == "__main__":
    main()
