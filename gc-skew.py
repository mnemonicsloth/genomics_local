import os, sys

genomeDir = "/home/jd/src/file-of-genomes/done/fna.d/"

fadir = '/home/jd/src/file-of-genomes/species.fa/'

chunk_size = 10000


def read_fna_file(fna_file):
    buf = ""
    for line in open(os.path.join(genomeDir, fna_file), 'r'):
        if line.startswith(">"):
            pass
        else:
            buf += line.rstrip("\n")
    return buf


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

    
# def cumulativize(tally):
#     acc = []
#     sum = (0,0)
#     for gc in tally:
#         sum = (sum[0] + gc[0], sum[1] + gc[1])
#         acc.append(sum)
#     return acc


def accumulate(list):
    acc = []
    sum = 0
    for n in list:
        sum += n
        acc.append(sum)
    return acc


def gcskew(g, c):
    return (g-c)/(g+c)
    

def read_fa_file(filename):
    acc = []
    for line in open(os.path.join(fadir, species, filename), 'r'):
        if line.startswith(">"):
            bits = line.strip("> +-\n").split("_")
            gene = "_".join(bits[:2])
            start = bits[-2]
            end = bits[-1]
            acc.append(  (gene, int(start), int(end))  )
        else:
            pass
    return acc


def partition_fa_list(fa_list, ori, term):
    end = lambda x: x[2]
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
    pos = gene[2]
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




species = sys.argv[1]


files = os.listdir(os.path.join(fadir, species))

report = []
i = 0
for f in files:
    fna_file = f.split(".fa")[0] + ".fna"
    buf = read_fna_file(fna_file)
    tally = [ do_chunk(c) for c in partition_chunks(buf, chunk_size)]
    local_gc_skew = [ gcskew(g,c) for g,c in tally ]
    cgc_skew = accumulate(local_gc_skew)
    ori = 10000 * cgc_skew.index(min(cgc_skew))
    term = 10000 * cgc_skew.index(max(cgc_skew))
    n = len(buf)
    report.append((f, ori, term, n))





    fa_list = read_fa_file(f)
    l,m,r = partition_fa_list(fa_list, ori, term)
    Lrep, Rrep = make_replichores(l,m,r, ori, term)
    data = []
    for gene in Lrep:
        dist = rdist(gene, ori, term, n)
        length = gene[2] - gene[1]
        data.append((dist, length, 0))
    for gene in Rrep:
        dist=rdist(gene, ori, term, n)
        length = gene[2] - gene[1]
        data.append((dist, length, 1))
    data.sort()
    fields = ['origin_dist', 'gene_length', 'replichore']
    filename = 'gene-dist.csv'
    with open(filename, 'w') as csvfile:
        w = csv.writer(csvfile)
        w.writerow(fields)
        w.writerows(data)
        w.close()
    













