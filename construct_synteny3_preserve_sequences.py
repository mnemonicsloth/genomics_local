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



def buildSyntenicReplichores(dir):
    # SR: syntenic replichore representation
    sr = {}
    d = process_families_file(dir)
    for file, genome in d.items():
        oriterm = getOriTerm(file, d)
        sr[file] = {}
        # print(oriterm)
        # print(list(d[file].keys()))
        # print('\n\n')
        for chromosomeID, geneSeq in genome.items():
            if chromosomeID in oriterm:
                ori, term, n = oriterm[chromosomeID]
                l,m,r = partition_fa_list(geneSeq, ori, term)
                Lrep, Rrep = make_replichores(l,m,r, ori, term)
                sr[file][chromosomeID] = (Lrep, Rrep)
    return d, sr


def get_chromosome_length(replichores):
    for geneSeq in replichores:
        largest = 0
        if len(geneSeq) == 0:
            continue
        n = max(geneSeq, key=lambda x: x[1][1]) [1][1]
        if n > largest:
            largest = n
    return largest





















# def getShortChromosome(sr):
#     shorties = {}
#     for filename, genome in sr.items():
#         min_n = 0
#         shortest = ''
#         for chromosomeID, replichores in genome.items():
#             n = get_chromosome_length(replichores)
#             if min_n == 0:
#                 min_n = n
#                 shortest = chromosomeID
#             else:
#                 if n < min_n:
#                     min_n = n
#                     shortest = chromosomeID
#         shorties[filename] = genome[shortest]
#     return shorties

        
# def getLongChromosome(sr):
#     longOnes = {}
#     for filename, genome in sr.items():
#         max_n = 0
#         longest = ''
#         for chromosomeID, replichores in genome.items():
#             n = get_chromosome_length(replichores)
#             if n > max_n:
#                 max_n = n
#                 longest = chromosomeID
#         longOnes[filename] = genome[longest]
#     return longOnes
    









# ====================================================
#
#            chromosome matching
#
# ====================================================



def listExtraChromosomes(sr):
    for k,v in sr.items():
        #if len(v) > 1:
        print(str(len(v)) + " " + k)


def getReplichoreLengths(replichorePair):
    L = len(replichorePair[0])
    R = len(replichorePair[1])
    return (L+R, L, R)

def getGenomeDimensions(genome):
    acc = []
    for replichorePair in genome.values():
        acc.append(getReplichoreLengths(replichorePair))
    return acc

def getSpeciesDimensions(sr):
    acc = []
    for filename, genome in sr.items():
        dim = getGenomeDimensions(genome)
        acc.append((filename, dim))
    return acc



def joinFileNames(name1, name2):
    return name1 + "$" + name2

def splitFileName(name):
    return tuple(name.split("$"))


def getChromosomeBySize(genome, size):
    for replichores in genome.values():
        if size == getReplichoreLengths(replichores)[0]:
            return replichores
    return None


def getBiggestChromosomes(sr):
    replichorePairs = {}
    chromosomeIDs = {}
    for filename, genome in sr.items():
        maxChromoName, maxReplichorePair = max(genome.items(),
                                             key=lambda x: getReplichoreLengths(x[1]))
        replichorePairs[filename] = maxReplichorePair
        chromosomeIDs[filename] = maxChromoName
    return replichorePairs, chromosomeIDs
        

def consolidateReplichores(replichorePair):
    L = {geneMarker for (geneMarker, (start,end)) in replichorePair[0]}
    R = {geneMarker for (geneMarker, (start,end)) in replichorePair[1]}
    genomeSet = L.union(R)
    return genomeSet


def calculateChromosomeOverlap(chr1, chr2):
    set1 = consolidateReplichores(chr1)
    set2 = consolidateReplichores(chr2)
    numerator = len(set1.intersection(set2))
    return min( numerator/len(set1), numerator/len(set2) )


def calculateOverlapsUsingReference(filename2pair):
    refpair = nthval(filename2pair, 0)
    for filename, pair in filename2pair.items():
        overlap = calculateChromosomeOverlap(refpair, pair)
        print(f"{filename}:     {overlap}")

        







# def consolidateReplichores(reps):
#     genomeSet = {}
#     for chromosomeID, geneseq in reps.items():
#         L = {geneMarker for (geneMarker, (start,end)) in geneseq[0]}
#         R = {geneMarker for (geneMarker, (start,end)) in geneseq[1]}
#         genomeSet[chromosomeID] = L.union(R)
#     return genomeSet






def matchChromosomes(ref, test):
    refset = consolidateReplichores(ref)
    testset = consolidateReplichores(test)

    if len(refset) > len(testset):
        reversedPairs = matchChromosomes(test, ref)
        return [(y,x) for x,y in reversedPairs]

    pairs = []
    for key, geneset in refset.items():
        bestkey = ''
        bestlen = 0
        for key2, geneset2 in testset.items():
            matchLength = len(geneset.intersection(geneset2))
            if matchLength > bestlen:
                bestkey = key2
                bestlen = matchLength
        if bestkey not in testset:
            raise Exception("matchChromosome: key2 not in testset.  No matching chromosome found")
        del testset[bestkey]
        pairs.append((key, bestkey))
    return pairs
        

def testChromosomeMatcher(genome1, genome2):
    f = lambda x: len(x[0])+len(x[1])
    return [(f(genome1[x]), f(genome2[y])) for x,y in matchChromosomes(genome1, genome2)]



# will give bad data if the reference genome has more chromosomes
# than a typical member of the set
# I'm not fixing this because we're already assuming the reference
# genome is a good choice
def matchAllChromosomes(sr):
    ref = nthval(sr,0)
    matches = {}
    sequence = [x for x,y in matchChromosomes(ref, ref)]

    assoc = lambda x, xs: [b for a,b in xs if x==a][0]
    
    for file,genome in sr.items():
        pairs = matchChromosomes(ref, genome)

        # if we have a multiple chromosome failure this is where it will happen
        orderedChromosomeLabels = [assoc(x, pairs) for x in sequence]
        matches[file] = orderedChromosomeLabels
    return matches


def buildChromosomeArray(sr):
    ref = nthkey(sr, 0)
    matches = matchAllChromosomes(sr)
    chromosomeArray = []
    for i in range(len(matches[ref])):
        chromosomeSequences = {}
        for file, table in matches.items():
            chromosomeSequences[file] = sr[file][table[i]]
        chromosomeArray.append(chromosomeSequences)
    return chromosomeArray





# ====================================================
#
#            replichore matching
#
# ====================================================



def matchReplichores(ref, test):
    refL = { geneMarker for (geneMarker, (start, end)) in ref[0]}
    refR = { geneMarker for (geneMarker, (start, end)) in ref[1]}

    testL = { geneMarker for (geneMarker, (start, end)) in test[0]}
    testR = { geneMarker for (geneMarker, (start, end)) in test[1]}

    LL = refL.intersection(testL)
    LR = refL.intersection(testR)
    RL = refR.intersection(testL)
    RR = refR.intersection(testR)

    return list(map(len, [LL, LR, RL, RR]))


def scoreReplichoreMatch(quad):
    n = sum(quad)
    return (quad[0] + quad[3])/n


def pairOffReplichores(ref, test):
    score = scoreReplichoreMatch(matchReplichores(ref, test))
    if score >= .5:
        return (test[0], test[1])
    else:
        return (test[1], test[0])
    

def compareReplichores(replichores, ref):
    acc = []
    for file, seq in replichores.items():
        match = matchReplichores(replichores[ref], seq)
        acc.append(  (file, match, scoreReplichoreMatch(match))  )
    acc.sort(key=lambda x: x[2])
    return acc


def scoreReplichoreReference(replichores, ref):
    acc = compareReplichores(replichores, ref)
    goodMatches = [(file, match, score) for file, match, score in acc
                   if score > .75 or score < .25]
    return len(goodMatches) / len(acc)


def rateReferences(replichores):
    keys = list(replichores.keys())
    scores = [(k, scoreReplichoreReference(replichores, k)) for k in keys]
    scores.sort(key= lambda x: x[1])
    scores.reverse()
    return scores


def findBestReference(replichores):
    return rateReferences(replichores)[0][0]


def splitChromosomeIntoReplichores(dictOfChromosomes):
    reffilename = findBestReference(dictOfChromosomes)
    ref = dictOfChromosomes[reffilename]
    left = {}
    right = {}
    for filename, replichores in dictOfChromosomes.items():
        match = pairOffReplichores(ref, replichores)
        left[filename] = match[0]
        right[filename] = match[1]
    return left, right, reffilename



def buildReplichoreArray(chromosomeArray):
    ra = []
    for dictOfChromosomes in chromosomeArray:
        ra.append(splitChromosomeIntoReplichores(dictOfChromosomes))
    return ra

























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

    aboth = [a for a in A if a in Both]
    bboth = [b for b in B if b in Both]
    n = len(Both)
    buf = [-1]*n

    for i in range(n):
        j = bboth.index(aboth[i])
        buf[j] = i
    if -1 in buf:
        raise Exception(buf)
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
        return []

    safety = [x for x in gnl if type(x) == type("foo")]
    if safety != []:
        raise Exception(safety)

    n = len(gnl)
    
    def getdir(curr, prev):
        diff = curr - prev
        if diff in {1, -1}:
            return diff
        elif curr == 0 and prev == n-1:
            return 1
        elif curr == n-1 and prev == 0:
            return -1
        else:
            return 0

    count = 1
    seqs = []
    for i in range(1,n):
        dir = getdir(gnl[i], gnl[i-1])
        if dir in {1, -1}:
            count += 1
        else:
            seqs.append(count)
            start = i
            count = 1
    seqs.append(count)

    if getdir(gnl[0], gnl[-1]) in {1,-1} and len(seqs) > 1:
        seqs[0] += seqs[-1]
        return seqs[:-1]
    else:
        return seqs


def backbone_stability(gnl):
    n = len(gnl)
    if n == 0:
        return 0
    def prev(i):
        return (i-1) % n
    def next(i):
        return (i+1) % n
    def neighbors(i):
        return (prev(i), next(i))

    # print([neighbors(gnl[i]) for i in range(n)])

    backbone_links = 0

    for i in range(n):
        if gnl[ prev(i) ] in neighbors(gnl[i]):
            backbone_links += 1

    return backbone_links / n


def process_directory(dir):
    d = process_families_file(dir)
    ref = get_reference_genome(dir)

    hist = []
    for k, v in d.items():
        matches = match_chromosomes(d[ref], v)
        syntenic_units = []
        bs_tally = []
        for c1, c2 in matches:
            diffs = compare_chromosomes(d[ref][c1], d[k][c2])
            segs = get_syntenic_segments(diffs)
            bs = backbone_stability(diffs)
            # print(f'{syntenic_units},   {segs}')
            syntenic_units += segs
            bs_tally.append((len(segs), bs))
        bs = sum([n*bs for n, bs in bs_tally]) / sum([n for n,bs in bs_tally])
        hist.append((syntenic_units, bs))
    return hist


# def compile_synteny_by_species():
#     dirs = sorted(os.listdir('ccout'))
#     output = []
#     for d in dirs:
#         d_full_path = os.path.join('ccout', d)
#         output.append((d, process_directory(d_full_path)))
#     return output


    
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

def count_singletons(xs):
    return sum([1 for x in xs if x==1])


def display(results):
    for species, counts in results.items():
        n = len(counts)
        
        synteny_lengths = [l for lengths, bs in counts for l in lengths]

        bs_acc = [bs for lengths, bs in counts]
        n_syntenic_units = [len(s) for s,bs in counts]
        
        n_mu = statistics.mean(n_syntenic_units)
        n_sigma = statistics.stdev(n_syntenic_units)

        bs_mu = statistics.mean(bs_acc)
        bs_sigma = statistics.stdev(bs_acc)
            
        l_mu = statistics.mean(synteny_lengths)
        l_sigma = statistics.stdev(synteny_lengths)

        singles = [count_singletons(lengths) for lengths, bs in counts]
        singles_mu = statistics.mean(singles)

        # print(f'{species:40} {float(singles_mu): 7.4} {float(n_mu):7.4} {float(bs_mu):7.4} {float(l_mu):7.4} {n:11}')

        print(f'{species}\t{float(n_mu):.4}\t{float(bs_mu):.4}\t{n}')

        
