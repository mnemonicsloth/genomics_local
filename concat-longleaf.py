#! /home/jd/anaconda3/bin/python

import os
import sys

# align_dir = sys.argv[1]
# fa_dir = sys.argv[2]

similarity_threshold = 0.95

output_name = sys.argv[1]

align_dir = os.path.join("/nas/longleaf/home/jdeszyck/pine/align", output_name)
fa_dir =  os.path.join("/nas/longleaf/home/jdeszyck/pine/species.fa", output_name)
output_dir =  "/nas/longleaf/home/jdeszyck/pine/concatenates"



def pick_reference_genome(directory):
    return os.listdir(directory)[0]

def get_ref_gene_seq(fn):
    ref_path = os.path.join(fa_dir, fn)
    return [name.strip("> +-\n") for name in open(ref_path, 'r') if name.startswith('>')]

def filename_fa_to_prot(fn):
    return fn.rstrip('fa') + 'prot'

def get_ref_gene_names():
    ref = pick_reference_genome(fa_dir)
    names = get_ref_gene_seq(ref)
    return [filename_fa_to_prot(ref) + '&' + name for name in names]


ref_seq = get_ref_gene_names()
ref_set = set(ref_seq)

fam_files = [f for f in os.listdir(align_dir) if f.endswith('.fa')]



def read_fam_file(filename):
    fn = os.path.join(align_dir, filename)
    genes = {}
    gene_name = ''
    gene_text = ''
    for line in open(fn, 'r'):
        if line.startswith('>'):
            if gene_text == '' and gene_name != '':
                raise Exception('read_fam_file: found empty gene')
            elif gene_name != '':
                genes[gene_name] = gene_text
            gene_name = line.strip('>\n')
            gene_text = ''
        else:
            gene_text += line.strip()
    genes[gene_name] = gene_text
    return genes



filenames = [filename_fa_to_prot(f) for f in os.listdir(fa_dir)]
genes = {}
analogs = {'nil': []}


def process_fam_file_data(dict):
    for k,v in dict.items():
        genes[k] = v

    nonref = [name for name in dict.keys() if name not in ref_set]
    ref    = [name for name in dict.keys() if name in ref_set]

    if len(ref) == 0:
        analogs['nil'].append(nonref)
    elif len(ref) > 1:
        raise Exception('process_fam_file_data: Too many reference genes in fam file')
    else:
        ref = ref[0]
        analogs[ref] = nonref

        
for f in fam_files:
    process_fam_file_data(read_fam_file(f))



ref_seq = [gene for gene in ref_seq if gene in analogs]


def getfn(gene):
    s =  gene.split("&")[0]
    return s

def concatenate():
    c = {}
    for f in filenames:
        c[f] = ""
    i = 0
    index = []

    for r in ref_seq:
        # fn = getfn(r)
        # c[fn] = c[fn] + genes[r]
        n = len(genes[r])
        gene_start = i
        gene_end = i + n - 1

        c[filenames[0]] = c[filenames[0]] + genes[r]
        
        for f in filenames[1:]:
            anlg = [marker for marker in analogs[r] if marker.startswith(f)]
            if len(anlg) == 0:
                c[f] = c[f] + n*'-'
            elif len(anlg) > 1:
                raise Exception("concatenate: specimen occurs more than once")
            else:
                c[f] = c[f] + genes[anlg[0]]

        i = i+n
        index.append((r, gene_start, gene_end))
    return c, index



def output_index(idx):
    output_file = os.path.join(output_dir, output_name + ".index.txt")
    with open(output_file, "w") as f:
        for line in idx:
            f.write(" ".join(str(x) for x in line))
            f.write("\n")

def output_gene_equivalences():
    gene_equivalences = [[r] + analogs[r] for r in ref_seq]
    output_file = os.path.join(output_dir, output_name + ".gene_equivalences.txt")
    with open(output_file, "w") as f:
        for line in gene_equivalences:
            f.write(" ".join(x for x in gene_equivalences[0]))


def remove_divergent_genes(conc):
    ref = conc[filenames[0]]
    for f in filenames[1:]:
        seq = conc[f]
        n = len(ref)
        if n != len(seq):
            raise Exception("remove_divergent_genes: Error: genomes of different lengths")
        diffs = 0
        for i in range(n):
            if ref[i] != seq[i]:
                diffs = diffs + 1
        pct_sim = (n-diffs)/n
        if pct_sim < similarity_threshold:
            del conc[f]
    return conc
            


def output_fa(conc):
    output_file = os.path.join(output_dir, output_name + ".fa")
    with open(output_file, "w") as f:
        for file in filenames:
            if file in conc:
                # print("MATCHED: " + file)
                f.write(">" + file + "\n")
                f.write(conc[file] + "\n")
            # else:               
              #   print("SKIPPED: " + file)


def go():
    c, index = concatenate()
    output_index(index)
    output_gene_equivalences()
    output_fa(remove_divergent_genes(c))


go()






