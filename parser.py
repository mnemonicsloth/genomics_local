#! /usr/bin/python

import os
import sys
from os import listdir
from os.path import isfile, join, basename
import re


def rev_comp(seq):
        new_seq=""
        SEQ = seq[::-1]
        SEQ = SEQ.upper()
        complement = {"C":"G", "G":"C", "A":"T", "T":"A",\
                      "R":"Y", "Y":"R", "K":"M", "M":"K",\
                      "S":"S", "W":"W", "N":"N", "B":"V",\
                      "V":"B", "D":"H", "H":"D"}
        for B in SEQ:
                if B not in complement:
                        print("!!!", B)
                else:
                        new_seq += complement[B]
        return new_seq


def translate(seq):
        d={}
        d["TTT"],d["TTC"],d["TTA"],d["TTG"]="F","F","L","L"
        d["CTT"],d["CTC"],d["CTA"],d["CTG"]="L","L","L","L"
        d["ATT"],d["ATC"],d["ATA"],d["ATG"]="I","I","I","M"
        d["GTT"],d["GTC"],d["GTA"],d["GTG"]="V","V","V","V"
        d["TCT"],d["TCC"],d["TCA"],d["TCG"]="S","S","S","S"
        d["CCT"],d["CCC"],d["CCA"],d["CCG"]="P","P","P","P"
        d["ACT"],d["ACC"],d["ACA"],d["ACG"]="T","T","T","T"
        d["GCT"],d["GCC"],d["GCA"],d["GCG"]="A","A","A","A"
        d["TAT"],d["TAC"],d["TAA"],d["TAG"]="Y","Y","*","*"
        d["CAT"],d["CAC"],d["CAA"],d["CAG"]="H","H","Q","Q"
        d["AAT"],d["AAC"],d["AAA"],d["AAG"]="N","N","K","K"
        d["GAT"],d["GAC"],d["GAA"],d["GAG"]="D","D","E","E"
        d["TGT"],d["TGC"],d["TGA"],d["TGG"]="C","C","*","W"
        d["CGT"],d["CGC"],d["CGA"],d["CGG"]="R","R","R","R"
        d["AGT"],d["AGC"],d["AGA"],d["AGG"]="S","S","R","R"
        d["GGT"],d["GGC"],d["GGA"],d["GGG"]="G","G","G","G"
        i=0
        tmp=[]
        while i in range(len(seq)-3):
                codon= seq[i:i+3]
                if codon not in d.keys():
                        tmp.append("X")
                else:
                        tmp.append(d[codon])
                        if d[codon]=="*":
                                break
                i += 3
        prot="".join(tmp)
        return  prot

#note: silently drops the last letters if len(seq) is not a multiple of 3



# def my_translate(seq):
#         d={}
#         d["TTT"],d["TTC"],d["TTA"],d["TTG"]="F","F","L","L"
#         d["CTT"],d["CTC"],d["CTA"],d["CTG"]="L","L","L","L"
#         d["ATT"],d["ATC"],d["ATA"],d["ATG"]="I","I","I","M"
#         d["GTT"],d["GTC"],d["GTA"],d["GTG"]="V","V","V","V"
#         d["TCT"],d["TCC"],d["TCA"],d["TCG"]="S","S","S","S"
#         d["CCT"],d["CCC"],d["CCA"],d["CCG"]="P","P","P","P"
#         d["ACT"],d["ACC"],d["ACA"],d["ACG"]="T","T","T","T"
#         d["GCT"],d["GCC"],d["GCA"],d["GCG"]="A","A","A","A"
#         d["TAT"],d["TAC"],d["TAA"],d["TAG"]="Y","Y","*","*"
#         d["CAT"],d["CAC"],d["CAA"],d["CAG"]="H","H","Q","Q"
#         d["AAT"],d["AAC"],d["AAA"],d["AAG"]="N","N","K","K"
#         d["GAT"],d["GAC"],d["GAA"],d["GAG"]="D","D","E","E"
#         d["TGT"],d["TGC"],d["TGA"],d["TGG"]="C","C","*","W"
#         d["CGT"],d["CGC"],d["CGA"],d["CGG"]="R","R","R","R"
#         d["AGT"],d["AGC"],d["AGA"],d["AGG"]="S","S","R","R"
#         d["GGT"],d["GGC"],d["GGA"],d["GGG"]="G","G","G","G"
#         i=0
#         amino_acids=[]
#         while i + 2 in range(len(seq)):
#                 codon = seq[i:i+3]
#                 if codon not in d:
#                         amino_acids.append("X")
#                 else:
#                         amino_acids.append(d[codon])
#                         if d[codon] == "*":
#                                 break
#                 i += 3
#         protein = "".join(amino_acids)
#         return protein
                









fnadir = "fna.d"
gffdir = "gff.d"
resultsdir = "results"

fnas = [f for f in listdir(fnadir) if isfile(join(fnadir, f))]
gffs = [f for f in listdir(gffdir) if isfile(join(gffdir, f))]
fnas.sort()
gffs.sort()

fnaset = { f.rstrip(".fna") for f in fnas}
gffset = { g.rstrip(".gff") for g in gffs}

if fnaset != gffset:
        raise ValueError




def check_matching_stems(fnas, gffs):
        odd_ones = []
        for f,g in zip(fnas, gffs):
                if f.rstrip(".fna") != g.rstrip(".gff"):
                        odd_ones.append((g,f))





def process_gff(file):
        f = open(gffdir + "/" + file, "r")
        contents = []
        for l in f:
                if l[0] != "#":
                        entries=l.strip("\n").split("\t")
                        if entries[2] == "CDS":
                                resu = [entries[0], entries[3], entries[4], entries[6]]
                                contents.append(resu)
        f.close()
        return contents


def process_fna(fnafile):
        data = {}
        f = open(fnadir + "/" + fnafile, "r")
        for l in f:
                # have checked the data.  all have ">" on line 1.  so for now just
                # assume there is a starting contig
                if l[0] == ">":
                        header = l.strip(">").strip("\n").strip("\r").split(" ")
                        contig = header[0]
                        data[contig] = []
                else:
                        data[contig].append(l.strip("\n").strip("\r").upper())
        f.close()
        for contig in data:
                data[contig] = "".join(data[contig])
        return data



# Pass filename as a parameter so we can use it in
# an error message if needed.  We also need it to determine what
# filenames to output to.
def process_genes(file_stem, gffdata, fnadata):
    if len(gffdata) > 4:
        h = open(resultsdir + "/" + file_stem + ".fa", "w")
        g = open(resultsdir + "/" + file_stem + ".prot", "w")
        for entry in gffdata:
            contig = entry[0]
            start,end = int(entry[1]), int(entry[2])
            sense = entry[3]
            if start > end:
                print(file_stem + ": starting index > final index in contig: ", contig)
            if sense == "+":
                gene = fnadata[contig][start-1 : end]
            else:
                gene = rev_comp(fnadata[contig][start-1 : end])
            h.write(">" + contig + "_" + str(start) + "_" + str(end) + \
                    " " + sense + "\n" + gene + "\n")
            prot = translate(gene)
            if len(gene) - 3 == len(prot) * 3:
                g.write(">" + contig + "_" + str(start) + "_" + str(end) + \
                        " " + sense + "\n" + prot + "\n")
            # else:
              #  print(file_stem +":" + contig +  ": length of gene:" + str(len(gene) - 3) + ".  3x length of protein: " \
              #        + str(3 * len(prot)))
        h.close()
        g.close()
    else:
            print(file_stem,": gffdata too short, only ", len(gffdata), " entries.")
            




                


def main():
        for g,f in zip(gffs, fnas):
                stem = g.rstrip(".gff")
                gdata = process_gff(g)
                fdata = process_fna(f)
                process_genes(stem, gdata, fdata)

