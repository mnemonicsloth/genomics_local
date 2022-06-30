import os
import sys


arguments = sys.argv

if "&" in arguments:
	arguments.remove("&")

gene_path= sys.argv[-2]
if gene_path[-1]=="/":
	pass
else:
	gene_path += "/"
	
genome_path = sys.argv[-1]
if genome_path == "/":
	pass
else:
	genome_path += "/"


species=["SP"]

tmp=os.listdir(genome_path)

liste={}
strains={}
for sp in species:
	liste[sp]=[]
	strains[sp]=[]
	for stuff in tmp:
		if stuff.endswith(".prot"):
			strains[sp].append(stuff.split(".prot")[0])
			liste[sp].append(stuff)

#liste={}
#for sp in species:
#	liste[sp]=[]
#	for st in strains[sp]:
#		truc = st + '.prot'
#		liste[sp].append(truc)



print('Load IDs')

parent={}
seq={}
for sp in species:
	parent[sp]={}
	for st in strains[sp]:
		truc = st + '.fa'
		f=open(genome_path + truc ,'r')
		for l in f:
			if l[0] == '>':
				id = l.strip('>').rstrip('\r').split(' ')[0]
				parent[sp][id]=st
		f.close()





print 'Load proteins'

shortcut={}
robert={}
corres={}
for sp in species:
	corres[sp],robert[sp]={},{}
	fichiers = os.listdir(gene_path )
	for truc in fichiers:
		if truc.endswith('.prot.align'):
			ortho = truc.split('.prot')[0]
			truc2 = ortho + '.fa'
			if truc2 not in fichiers:
				robert[sp][ortho]={}
				corres[sp][ortho]={}
				print(truc)
				f=open(gene_path  + truc ,'r')
				for l in f:
					if l[0] == '>':
						id=l.strip('\n').lstrip('>').rstrip('\r')
						robert[sp][ortho][id]=''
						if 1==1:
							shortcut[id] = 'y'
						corres[sp][ortho][id]=id
					else:
						robert[sp][ortho][id]+=l.strip('\n').rstrip('\r')
				f.close()




print 'Load nucleotides'


tmp={}
for sp in species:
	fichiers = os.listdir(genome_path)
	for truc in fichiers:
		if truc.endswith('.fa'):
			f=open(genome_path + truc ,"r")
			for l in f:
				if l[0]==">":
					id = l.strip("\n").strip(">").rstrip('\r').split(' ')[0]
					#print '1 ',id
					id = truc.split(".fa")[0] + ".prot" + "&" + id
					tmp[id] = []
				elif shortcut.has_key(id):
					tmp[id].append( l.strip("\n").rstrip('\r').upper())
			f.close()



seq={}
for id in tmp:
	seq[id] = ''.join(tmp[id])


tmp=[]


print 'Back translate'




back={}
for sp in species:
	back[sp]={}
	for ortho in robert[sp]:
		#print 'Back translate ',sp,' ',ortho
		back[sp][ortho]={}
		if 1==1	:
			tag=0
			for id in robert[sp][ortho]:
				back[sp][ortho][id] = ""
				tag+=1
				resu = ""
				align = robert[sp][ortho][id]
				ID = corres[sp][ortho][id]
				i,j=0,0
				while i < len(align):
					aa = align[i]
					if aa == "-":
						resu += "---"
					else:
						codon = seq[ID][j:j+3]
						resu += codon
						j+=3
					if len(resu) == 60:
						back[sp][ortho][id] += resu
						resu=""
					i+=1
				if len(resu) >= 1:
					back[sp][ortho][id] += resu



for sp in species:
	print sp, ' ',back[sp].keys()



print "Write fasta"


for sp in species:
	print "write", sp
	for ortho in back[sp]:
		h=open(gene_path +  ortho + ".fa","w")
		for id in back[sp][ortho]:
			if 1==1:
				h.write(">" + id + "\n")
				i=0
				while i < len(back[sp][ortho][id]):
					h.write(back[sp][ortho][id][i:i+60] + "\n")
					i+=60
		h.close()


















