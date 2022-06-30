# for reference: synteny is a sequence of core genes in the order in
# which they appear in a genome.  This program takes the output of
# core cruncher and constructs a synteny (an ordered list of core
# genes ['famXXX', 'famYYY', 'famZZZ', ...]) for each genome in
# Yersinia_pestis/species.prot


import os, sys, math
from itertools import permutations
import statistics

working_dir = "ccout/Yersinia_pestis"

genome_files_dir = 'species.prot'
famfilename = 'families_core.txt'
output_dir = 'synteny.out'
f = open(os.path.join(working_dir,famfilename), 'r')

# returns
# ('GCF_000009065.1_ASM906v1_genomic.prot', 'fam2619', ('NC_003143.1', 3171247, 3172095))
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


def dist(x,y):
    return abs(x-y)


def dissimilarity(g1, g2, g1keyseq, g2keyseq):
    d = 0
    # print(f'{g1keyseq}  |||   {g2keyseq}')
    for i in range(len(g1keyseq)):
        d += abs(  len(g1[g1keyseq[i]]) - len(g2[g2keyseq[i]])  )
    return d


def match_chromosomes(g1, g2):
    # a genome is a dict whose keys are chromosome IDs and whose
    # values are ordered lists of items ('fam111', (start end)).
    #
    # returns a list of (chromosome_id1, chromosome_id2)
    
    k1 = sorted(list(g1.keys()), key=lambda x: len(g1[x]))
    k2 = sorted(list(g2.keys()), key=lambda x: len(g2[x]))

    # dissimilarity is only defined for seqs of the same length 
    # def dissimilarity(g1_key_seq, g2_key_seq):
    #     d = 0
    #     for i in range(len(g1_key_seq)):
    #         d += abs( len(g1[g1_key_seq[i]]) - len(g2[g2_key_seq[i]])  )
    #     return d

    if len(g1) == len(g2):
        return list(zip(k1, k2))

    if len(g1) < len(g2):
        matching = min(permutations(g2, len(g1)), key=lambda p: dissimilarity(g1, g2, k1, p))
        return list(zip(k1, matching))

    if len(g1) > len(g2):
        matching = min(permutations(g1, len(g2)), key=lambda p: dissimilarity(g1, g2, p, k2))
        return list(zip(matching, k2))


def list_chromosomes(genome):
    for k in genome.keys():
        print(k, len(genome[k]))

def list_genomes(d):
    for k in d.keys():
        print(k)
        for l in d[k].keys():
            print("     ", l, len(d[k][l]))
    
# def detect_overlaps(chromosome_list)     assuming none for now1


def extract_gene_sequence(chromosome_list):
    return [gene for gene, (start, end) in chromosome_list]


def do_chromo_sets(a, b):
# this takes geneID lists, not (geneID, (start, end)) lists
    A = set(a)
    B = set(b)
    return A-B, B-A, A.intersection(B)


def compare_chromosomes(ca, cb):
    A = extract_gene_sequence(ca)
    B = extract_gene_sequence(cb)
    Aonly, Bonly, Both = do_chromo_sets(A, B)

    buf = ['_'] * len(B)
    for i in range(len(A)):
        if A[i] in Aonly:
            pass
        else:
            j = B.index(A[i])
            buf[j] = i
    while '_' in buf:
        buf.remove('_')
    return buf


 # def get_short_chromosomes(d):
 #    keys = list(d.keys())
 #    c = [0] * len(keys)
 #    for i in range(len(keys)):
 #        genome = d[keys[i]]
 #        shortest_key = min(genome.keys(), key=lambda x: len(genome[x]))
 #        c[i] = genome[shortest_key]
 #    return c
















 
def get_syntenic_segments(gene_number_list):
    gnl = gene_number_list
    if gnl == []:
        return 0

    safety = [x for x in gnl if type(x) == type("foo")]
    if safety != []:
        raise Exception(safety)
    
    low = min(gnl)
    high = max(gnl)

    start = 0
    count = 0
    seqs = []

    def getdir(curr, prev):
        diff = curr - prev
        if diff in {1, -1}:
            return diff
        elif curr == low and prev == high:
            return 1
        elif curr == high and prev == low:
            return -1
        else:
            return 0

    for i in range(1,len(gnl)):
        dir = getdir(gnl[i], gnl[i-1])
        if dir in {1, -1}:
            count += 1
        else:
            seqs.append((start, start+count))
            start = i
            count = 0
    seqs.append((start, start + count))

    nseqs = len(seqs)
    if getdir(gnl[0], gnl[-1]) in {1,-1} and nseqs > 1:
        nseqs -= 1
    return nseqs


    









def process_directory(dir):
    d = process_families_file(dir)
    # so now we have dict[file][chromosome] -> sequence-of-IDs-on-chromosome
    ref = get_reference_genome(dir)

    hist = []
    for k, v in d.items():
        matches = match_chromosomes(d[ref], v)
        syntenic_units = 0
        for c1, c2 in matches:
            diffs = compare_chromosomes(d[ref][c1], d[k][c2])
            syntenic_units += get_syntenic_segments(diffs)
        hist.append(syntenic_units)
    return hist


# def compile_synteny_by_species():
#     dirs = sorted(os.listdir('ccout'))
#     output = []
#     for d in dirs:
#         d_full_path = os.path.join('ccout', d)
#         output.append((d, process_directory(d_full_path)))
#     return output


def get_core_genome_size(dir):
    d_full_path = os.path.join("ccout", dir, "core")
    core_genome = os.listdir(d_full_path)
    return len(core_genome)



def get_core_genome_sizes():
    sizes = {}
    for dir in os.listdir("ccout"):
        sizes[dir] = get_core_genome_size(dir)
    return sizes






def compile_synteny_by_species():
    dirs = sorted(os.listdir('ccout'))
    dirs.remove('Borreliella_burgdorferi')
    acc = {}
    for d in dirs:
        print(d)
        d_full_path = os.path.join('ccout', d)
        acc[d] = process_directory(d_full_path)
    return acc







# def compile_species_synteny(dir):
#     print(dir)
#     d_full_path = os.path.join('ccout', dir)
#     f = open(os.path.join('synteny.out', dir), 'w')
#     output = process_directory(d_full_path)
#     f.write(str(output))
#     f.close()
#     return True


def tabulate_results(acc):
    tally = []
    for k,v in acc.items():
        n = len(v)
        mean = sum(v)/n
        devs = [(x-mean)**2 for x in v]
        var = sum(devs)/(n-1)
        stdev = math.sqrt(var)
        tally.append( (k, mean, stdev, n) )
    return tally



def tabulate_normalized_results(acc, sizes):
    tally = []
    for k,v in acc.items():
        n = len(v)
        mean = sum(v)/n
        devs = [(x-mean)**2 for x in v]
        var = sum(devs)/(n-1)
        stdev = math.sqrt(var)

        size = sizes[k]
        mean = mean/size
        stdev = stdev/size
        tally.append( (k, mean, stdev, n) )
    return tally    
    






def display(results):
    for species, mean, stdev, n in results:
        print(f'{species:40} {float(mean):7.4} {stdev:7.4} {n:7}')



def save_output(results):
    outfile = open("normalized-synteny-scores.txt", "w")
    outfile.write("species\tmean\tstdev\tn\n")
    for species, mean, stdev, n in results:
        outfile.write(f'{species:40} {float(mean):7.4} {stdev:7.4} {n:7}')
