

print("USAGE: python master.py  -in input_folder   -out  output_folder   OPTIONAL: -freq  frequency_across_genomes  -ref ref_genome  -id unique/combined   -ext .fa/.fasta/.prot/.faa   -list path_to_genome_list   -score identify_score   -length sequence_length_conservation  -restart yes/no") 
print("Use -h to see the different options\n")

import os
import sys

from datetime import datetime
startTime = datetime.now()

path =  "/Users/lbobay/Desktop/dossiers/pan/test/database/"

out_path = "/Users/lbobay/Desktop/dossiers/pan/test/output/"

REF = "NA"

arguments =   sys.argv

if "-h" in arguments:
	print("\nUSAGE:")
	print("-in      input folder")
	print("-out     output folder")
	print("\nOPTIONS: ")
	print("-freq     Minimum frequency of the gene across genomes to be considered core (default= 90%, an ortholog is considered a core gene even if it is missing in 10% of the set of genomes)")
	print("-score    Identity score used by usearch to define orthologs in % (Default= 70)")
	print("-length   Minimum sequence length conservation used by to define orthologs (default= 80%)")
	print("-ref      Reference genome (default: first genome in folder will be used as reference). If you want to specify the reference genome to use, speficy the name of the file in the folder (e.g. -ref genome1.prot)")
	print("-id       Type of gene IDs in output files. Choose 'unique' if the same gene IDs are not found in different genomes or 'combined' to combine genome ID & gene ID (default= 'combined').  ")
	print("-ext	     File extensions .fa/.fasta/.prot/.faa (default: .fa )")
	print("-list     Path to a file containing the list of genomes to analyze (default: none, all the genomes in the folder will be analyzed by default)")
	print("-restart  Restart analysis from scratch yes or no (default= no). If yes is chosen, the program will erase the usearch output files and relaunch usearch")
	exit()

if "-in" in arguments:
	i = arguments.index("-in")
	path = arguments[i+1]
	if path[-1] == "/":
		pass
	else:
		path += "/"

if "-out" in arguments:
	i = arguments.index("-out")
	out_path = arguments[i+1]
	if out_path[-1] == "/":
		pass
	else:
		out_path += "/"

if "-freq" in arguments:
	i = arguments.index("-freq")
	freq = arguments[i+1]
else:
	freq="90"

if "-ref" in arguments:
	i = arguments.index("-ref")
	REF = arguments[i+1]

if "-ext" in arguments:
	i = arguments.index("-ext")
	ext = arguments[i+1]
else:
	ext = ".fa"
	
if "-list" in arguments:
	i = arguments.index("-list")
	FILE = arguments[i+1]
else:
	FILE = "none"

if "-score" in arguments:
	i = arguments.index("-score")
	score= arguments[i+1]
else:
	score="70"

if "-length" in arguments:
	i = arguments.index("-length")
	length = arguments[i+1]
else:
	length = "80"

if "-id" in arguments:
	i = arguments.index("-id")
	IDENTIFIANTS = arguments[i+1]
else:
	IDENTIFIANTS= "combined"

if "-restart" in arguments:
	i = arguments.index("-restart")
	restart = arguments[i+1]
else:
	restart= "no"

print("input = ", path)
print("output = ", path)



try:
	os.mkdir(out_path)
except OSError:
	pass

try:
	os.mkdir(out_path + "BBH")
except OSError:
	pass

if restart == "yes":
	os.system("rm " + out_path + "BBH/*" )	


try:
	os.mkdir(out_path + "core")
except OSError:
	os.system("rm -r " + out_path + "core" )	
	os.mkdir(out_path + "core")

tmp = os.listdir(path)


extensions=[".fa",".fasta",".prot",".faa"]
if ext in extensions:
	extensions.remove(ext)

memo=[]
tag=0
files = []
for stuff in tmp:
	if stuff.endswith(ext):
		files.append(stuff)
		tag=1
	else:
		for point in extensions:
			if stuff.endswith(point):
				if point not in memo:
					memo.append(point)

if len(memo) > 0:
	print("\n######################################")
	print("WARNING: other files have been detected with the following extensions: "," ".join(memo))
	print("The analysis will be conducted on the ",ext, " files")
	print("If you want to analyze other files, specify the file extensions to use with the option -ext  (e.g.  -ext .fasta)")
	print("######################################\n")

if tag==0:
	print("\n######################################")
	print("No files have been found with extension " + ext)
	print("Please specify the file extension with option -ext (e.g. -ext .prot)")
	print("Exiting...")
	print("######################################\n")
	exit()

toto,seq,size={},{},{}
check=[]
nb=0
for stuff in files:
	if 1==1:
		nb+=1
		if nb <= 3 or nb >= len(files) - 1:
			f=open(path + stuff ,"r")
			for l in f:
				if l[0]==">":
					id = l.strip(">").strip("\n").split(" ")[0]
					id = stuff + "&" + id
					toto[id]=[]
				else:
					toto[id].append(l.strip("\n"))
			f.close()
			for id in toto:
				seq[id] = "".join(toto[id]).upper()
				size[id] = len(seq[id])
				break
			SEQ = seq[id]
			A = SEQ.count("A")
			C = SEQ.count("C")
			G = SEQ.count("G")
			T = SEQ.count("T")
			dash = SEQ.count("-")
			tot = A + C + G + T + dash
			ratio = float(tot) / size[id]
			check.append(ratio)

#print("Check= ",check

if min(check) < 0.5 and max(check) < 0.5 :
	print("TYPE= Protein sequences")
	TYPE="protein"
elif min(check) > 0.5 and max(check) > 0.5 :
	print("TYPE= DNA sequences")
	TYPE="DNA"
elif min(check) < 0.5 and max(check) > 0.5 :
	print("PROBLEM: BOTH DNA and protein sequences were given. Exiting...")
	exit()

print("Building core genome with ",len(files)," genomes")

if REF =="NA":
	REF= files[0]

print("Reference genome= ", REF)

print(datetime.now() - startTime)


os.system("python  usearch_core.py  " + path + " " + out_path + " " + REF + " " + TYPE + " " + ext + " " + FILE + " " + freq + " " + score + " " + length + " " + IDENTIFIANTS)


print("Usearch complete")
print(datetime.now() - startTime)

os.system("python  big_cruncher.py  " + path + " " + out_path + " " + REF + " " + TYPE + " " + ext + " " + FILE + " " + freq + " " + score + " " + length + " " + IDENTIFIANTS)


print(datetime.now() - startTime)














