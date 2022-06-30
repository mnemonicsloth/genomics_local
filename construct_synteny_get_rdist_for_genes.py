# for reference: synteny is a sequence of core genes in the order in
# which they appear in a genome.  This program takes the output of
# core cruncher and constructs a synteny (an ordered list of core
# genes ['famXXX', 'famYYY', 'famZZZ', ...]) for each genome in
# Yersinia_pestis/species.prot


import os, sys
from itertools import permutations
import statistics

working_dir = "ccout/Yersinia_pestis"

genome_files_dir = 'species.prot'
famfilename = 'families_core.txt'
fnaDir = '/home/jd/src/file-of-genomes/done/fna.d/'





# f = open(os.path.join(working_dir,famfilename), 'r')







def prarray(arr):
    for x in arr:
        print(x)

def nthkey(d, n):
    return list(d.keys())[n]

def nthval(d, n):
    return list(d.values())[n]

def nthitem(d, n):
    return list(d.items())[n]








# returns
# ('GCF_000009065.1_ASM906v1_genomic.prot', 'fam2619', ('NC_003143.1', 3171247, 3172095)
def make_gene_marker(core_gene_id, gene_loc_str):
    file, right_side = tuple(gene_loc_str.split('&'))
    fragments = tuple(right_side.split('_'))
    end = fragments[-1] 
    start = fragments[-2]
    rest = fragments[0:-2]
    if len(rest) == 1:
        rest = rest[0]
    else:
        rest = '_'.join(rest)
    loc = (rest, int(start), int(end))
    return file, core_gene_id, loc

# modifies dict to include entries like
# ('GCF_000009065.1_ASM906v1_genomic.prot', [('fam2619', ('NC_003143.1', 3171247, 3172095))])
# net result: a dict of the form
#     filename ->  [ (geneid, loc), (geneid2, loc2), (geneid3, loc3)...]
def process_families_line(line, dict):
    line_entries = line.split('\t')
    markers = [make_gene_marker(line_entries[0], entry) for entry in line_entries[1:]]
    for file, gene, loc in markers:
        if file not in dict:
            dict[file] = []
        dict[file].append((gene, loc))

# returns a dict of the form
# chromosomeID -> coreGeneID, Startpos, Endpos
def assemble_by_chromosome(gene_list):
    chromosome_dict = {}
    for gene_id, (chromosome_id, start, end) in gene_list:
        if chromosome_id not in chromosome_dict:
            chromosome_dict[chromosome_id] = []
        chromosome_dict[chromosome_id].append(  (gene_id, (start, end))  )
    for chr_id, chr_gene_list in chromosome_dict.items():
        chromosome_dict[chr_id] = sorted(chr_gene_list, key=lambda x: x[1])
    return chromosome_dict


def process_families_file(dir):
    dict = {}
    for line in open(os.path.join(dir, famfilename), 'r'):
        process_families_line(line, dict)
    for file, data in dict.items():
        dict[file] = assemble_by_chromosome(data)
    return dict

# so now we have dict[file][chromosome] -> sequence-of-IDs-on-chromosome

def get_reference_genome(dir):
    return os.listdir(os.path.join(dir, genome_files_dir))[0]











def read_fna_file(fna_file):
    bufs = []
    chromosome = ''
    buf = ""
    for line in open(os.path.join(fnaDir, fna_file), 'r'):
        if line.startswith(">"):
            if chromosome != "":
                bufs.append((chromosome, buf))
                buf = ""
            chromosome = line.split(" ")[0].lstrip(">")
        else:
            buf += line.rstrip("\n")
    bufs.append((chromosome, buf))
    return bufs


def read_corresponding_fna(prot_file):
    fna_file = prot_file.split(".prot")[0] + ".fna"
    return read_fna_file(fna_file)



def do_chunk(chunk):
    g = 0
    c = 0
    for i in chunk:
        if i == "G" or i =="g":
            g += 1
        if i == "C" or i == "c":
            c += 1
        else:
            pass
    return (g, c)


def partition_chunks(buf, n):
    for i in range(0, len(buf), n):
        yield buf[i:i+n]


def accumulate(list):
    acc = []
    sum = 0
    for n in list:
        sum += n
        acc.append(sum)
    return acc


def gcskew(g, c):
    if g == c == 0:
        return 1
    else:
        return (g-c)/(g+c)



def partition_fa_list(fa_list, ori, term):
    end = lambda x: x[1][1]
    left = []
    middle = []
    right = []
    if ori < term:
        a=ori
        b=term
    else:
        a=term
        b=ori
    for i in fa_list:
        if end(i) < a:
            left.append(i)
        elif a <= end(i) <= b:
            middle.append(i)
        else:
            right.append(i)
    return left,middle,right



def rdist(gene, ori, term, n):
    # the Replichore DISTance function
    pos = gene[1][1]
    if ori < term:
        if pos < ori:
            return ori-pos
        elif ori <= pos <= term:
            return pos - ori
        else:
            return ori + (n-pos)
    else: # ie if term < ori
        if pos < term:
            return pos + n - ori
        elif term <= pos <= ori:
            return ori - pos
        else: # ie if ori < pos
            return pos - ori


def make_replichores(left,middle,right, ori, term):
    if ori < term:
        return left[::-1] + right[::-1], middle
    else:
        return middle[::-1], right + left



# filename must be a key in species_dict.  So it's a .prot file.
def getOriTerm(filename, species_dict):
    chunk_size = 10000
    bufs = read_corresponding_fna(filename)
    revised_bufs = [(head, text) for head, text in bufs
                    if head in species_dict[filename].keys()
                    # and len(text) > 100000
                    ]
    if revised_bufs == []:
        raise Exception(filename)
    oriTerms = {}
    for chromosomeID, text in revised_bufs:
        tally = [do_chunk(c) for c in partition_chunks(text, chunk_size)]
        local_gc_skew = [gcskew(g,c) for g,c in tally]
        cgc_skew = accumulate(local_gc_skew)
        ori = 10000 * cgc_skew.index(min(cgc_skew))
        term = 10000 * cgc_skew.index(max(cgc_skew))
        oriTerms[chromosomeID] = ( ori, term, len(text))
    return oriTerms





stems = os.listdir("ccout")

for stem in stems:
    # stem = "Escherichia_coli"    
    dir = os.path.join("ccout", stem)
    outfile_name = os.path.join("gene_locs", stem + "_loc.txt")
    outfile = open(outfile_name, "w")
    outfile.write("fam\tloc\n")
    outfile.flush()
    ref = get_reference_genome(dir)
    ffdict = process_families_file(dir)
    ot = getOriTerm(ref, ffdict)


    for chrid, (ori, term, n) in ot.items():
        refgenomedata = ffdict[ref][chrid]
        for gene in refgenomedata:
            fam = gene[0]
            loc = rdist(gene, ori, term, n)
            # print(f"{fam}\t{loc}\n")
            outfile.write(f"{fam}\t{loc}\n")
            # outfile.flush()
    print(stem)
    outfile.close()

