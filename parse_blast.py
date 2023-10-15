#!/usr/bin/python3

# parse_blast.py

# Hsin-Ying Chang <hyhazelchang@gmail.com>
# v1 2023/09/18

# Usage: python3 /home/xinchang/pyscript/pyscript_xin/parse_blast.py --in_dir=/scratch/xinchang/cyano17/cyano17.01/blast/ --out_dir=/scratch/xinchang/cyano17/cyano17.01/parse_blast/ --identity=1 --hsp=1

import argparse
import os
import glob

def main():
    parser = argparse.ArgumentParser(
            description=("Parse the blast results"),
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--in_dir",
                        type=str,
                        default=None,
                        help="xml files directory (out format 5).")
    parser.add_argument("--out_dir",
                        type=str,
                        default=None,
                        help="Output directory, please provide absolute path")
    parser.add_argument("--orthomcl_bpo",
                        type=int,
                        default=0)
    parser.add_argument("--identity",
                        type=int,
                        default=0)
    parser.add_argument("--hsp",
                        type=int,
                        default=0)
    
    # Defining variables from input
    args = parser.parse_args()
    in_dir = args.in_dir
    out_dir = args.out_dir
    orthomcl_bpo = args.orthomcl_bpo
    identity = args.identity
    hsp = args.hsp

    # Find the blast files
    xmlfiles = glob.glob(in_dir + "*.xml")

    # Check output directory
    if not os.path.exists(out_dir):
        os.system("mkdir -p " + out_dir)

    for xmlfile in xmlfiles:
        filename = os.path.basename(xmlfile).split(".")[0]
        # Parse out the blast results
        lines = [line.rstrip("\n").lstrip(" ") for line in open(xmlfile)]

        # Get hit index
        hit_index = Hit_index(lines)

        # Basic blast infomation
        blast = []
        for num in range(len(lines)):
            if "<Iteration_query-def>" in lines[num]:
                query = lines[num].replace("<Iteration_query-def>", "").replace("</Iteration_query-def>", "")              
        for num in range(1, len(hit_index)):
            blast.append([str(num), query])

        right = 1
        while right < len(hit_index):
            identical_length = 0
            match_length = 0
            Hsps = []
            for num in range(hit_index[right - 1], hit_index[right]):      
                if "<Hit_def>" in lines[num]:
                    Hit_def = lines[num].replace("<Hit_def>", "").replace("</Hit_def>", "")
                    blast[right-1].append(Hit_def)
                if "<Hit_accession>" in lines[num]:
                    Hit_acc = lines[num].replace("<Hit_accession>", "").replace("</Hit_accession>", "")
                    blast[right-1].append(Hit_acc)
                if identity: 
                    identical_length += Identity(lines, num)[0] 
                    match_length += Identity(lines, num)[1]
                if hsp:
                    Hsps_list = Hsp(lines, num, Hsps)
            if identity:
                percent_identity = str(round(identical_length / match_length * 100, 2))
                blast[right-1].append(percent_identity)
            if hsp:
                blast[right-1].append("".join(Hsps_list))
            right += 1

        # Print out the information
        ## concat the info
        for num in range(len(blast)):
            blast[num] = ";".join(blast[num])

        ## print out
        if orthomcl_bpo:
            out = open(out_dir + filename + ".bpo", "w")
        else:
            out = open(out_dir + filename + ".txt", "w")
        for info in blast:
            out.write(info + "\n")
        out.close()

def Hit_index(lines):
    count = 0
    hit_index = []
    while count < len(lines):
        if lines[count] == "<Hit>":
            hit_index.append(count)
        count += 1
    hit_index.append(len(lines))
    return hit_index

def Identity(lines, num):
    hsp_identity, hsp_match_length = 0, 0
    if "<Hsp_identity>" in lines[num]:
        hsp_identity = int(lines[num].replace("<Hsp_identity>", "").replace("</Hsp_identity>", ""))
    if "<Hsp_align-len>" in lines[num]: 
        hsp_match_length = int(lines[num].replace("<Hsp_align-len>", "").replace("</Hsp_align-len>", "")) 
    return hsp_identity, hsp_match_length

def Hsp(lines, num, Hsps):
    if "<Hsp_num>" in lines[num]:
        Hsp_num = lines[num].replace("<Hsp_num>", "").replace("</Hsp_num>", "")
        query_start = lines[num + 4].replace("<Hsp_query-from>", "").replace("</Hsp_query-from>", "")
        query_end = lines[num + 5].replace("<Hsp_query-to>", "").replace("</Hsp_query-to>", "")
        subject_start = lines[num + 6].replace("<Hsp_hit-from>", "").replace("</Hsp_hit-from>", "")
        subject_end = lines[num + 7].replace("<Hsp_hit-to>", "").replace("</Hsp_hit-to>", "")
        Hsp = Hsp_num + ":" + query_start + "-" + query_end + ":" + subject_start + "-" + subject_end + "."
        Hsps.append(Hsp)
    return Hsps

if __name__ == '__main__':
    main()
