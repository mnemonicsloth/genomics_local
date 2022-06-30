def rev_comp( seq ):
	new_seq=""
	SEQ = seq[::-1]
	SEQ = SEQ.upper()
	for B in SEQ:
		if B == "G":
			new_seq+="C"
		elif B == "C":
			new_seq+="G"
		elif B == "A":
			new_seq+="T"
		elif B == "T":
			new_seq+= "A"
		elif B == "R":
			new_seq+="Y"
		elif B == "S" or B == "W" or B == "N":
			new_seq += B
		elif B == "Y":
			new_seq+="R"
		elif B == "K":
			new_seq += "M"
		elif B == "M":
			new_seq+="K"
		elif B == "B":
			new_seq+="V"
		elif B == "V":
			new_seq += "B"
		elif B == "D":
			new_seq+="H"
		elif B =="H":
			new_seq+="D"
		else:
			print "!!! ",B
	return new_seq


#######################################################
# Translator

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
	return	prot


import os

species=[]
f=open('../species.txt','r')
for l in f:
	species.append(l.strip('\n'))


f.close()


#species = list(species[:100])



print species

PATH = '../results/'



liste={}
for sp in species:
	tmp= os.listdir(PATH + sp + '/genomes/')
	liste[sp]=[]
	nb=0
	for truc in tmp:
		if truc.endswith('.gff')  :
			nb+=1
			if nb >= 0: # and nb <= 5000:
				liste[sp].append(truc)
	print sp,' : ',nb,' genomes total'


dico={}
for sp in species:
	print 'Load GFF ',sp
	dico[sp]={}
	for file in liste[sp]:
		dico[sp][file.rstrip(".gff")] = []
		f=open(PATH + sp + '/genomes/' + file,"r")
		for l in f:
			if l[0] != "#":
				a=l.strip("\n").split("\t")
				if a[2] ==  "CDS":
					resu = [a[0],a[3],a[4],a[6]]
					dico[sp][file.rstrip(".gff")].append(resu)
		f.close()



#print dico[sp][file.rstrip(".gff")][:3]

print 'LOAD sequences'

tmp={}
for SP in species:
	print 'Load ', SP
	tmp[SP]={}
	for file in liste[SP]:
		sp = file.rstrip(".gff")
		#print sp
		tmp[SP][sp] = {}
		f=open(PATH + SP + '/genomes/' + sp + ".fna","r")
		for l in f:
			if l[0] == ">":
				a=l.strip(">").strip("\n").strip("\r").split(" ")
				contig =  a[0]
				tmp[SP][sp][contig] = []
			else:
				tmp[SP][sp][contig].append(l.strip("\n").strip("\r").upper())
		f.close()


seq={}
for SP in species:
	seq[SP]={}
	for sp in tmp[SP]:
		seq[SP][sp]={}
		for contig in tmp[SP][sp]:
			seq[SP][sp][contig] = "".join(tmp[SP][sp][contig])


tmp=''


#print seq[SP]['Gapicola_M6_G3']


redo = []

for SP in species:
	done = os.listdir(PATH + SP + '/genes/')
	NB=0
	for sp in dico[SP]:
		NB+=1
		if len(dico[SP][sp]) > 4 and sp + ".fa" not in done:
			print 'Write ', SP,' ',NB
			h=open(PATH + SP + '/genes/'  + sp + ".fa" , "w")
			g=open(PATH + SP + '/genes/' + sp + ".prot" , "w")
			for resu in dico[SP][sp]:
				contig = resu[0]
				deb,fin = int(resu[1]),int(resu[2])
				sens=resu[3]
				if deb > fin:
					print "problem ",sp," ",resu
				if sens== "+":
					gene = seq[SP][sp][contig][deb-1 : fin]
				else:
					gene = rev_comp(seq[SP][sp][contig][deb-1 : fin])
				h.write(">" + contig + "_" + str(deb) + "_" + str(fin) + " " + sens + "\n" + gene + "\n")
				prot = translate(gene)
				if len(gene) - 3 == len(prot) * 3:
					g.write(">" + contig + "_" + str(deb) + "_" + str(fin) + " " + sens + "\n" + prot + "\n")
				#else:
					#print len(gene) - 3," ",len(prot) * 3
			f.close()
			h.close()
			#g.close()
			gcf = sp + ".fa"
			#if robert[SP].has_key(gcf):
				#print 'Erase ', gcf
				#os.system('rm ../genes/' + robert[SP][gcf])
		else:
			print SP,' ',sp,' ',NB,' empty or already done'
















