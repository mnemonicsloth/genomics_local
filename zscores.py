import os
import numpy as np
import scipy.stats as ss

datadir = "sliding_window"



def destringify(line):
    tofloats = [1,2,3,4]
    for i in tofloats:
        line[i] = float(line[i])
    return line


def summarize(data):
    cleandata = [d for d in data if not d==float('nan')]
    zs = ss.zscore(cleandata, nan_policy="omit")
    zmin = min(zs)
    wheremin = np.where(zs==zmin)
    zmax = max(zs)
    wheremax = np.where(zs==zmax)
    avg = np.mean(cleandata)
    med = np.median(cleandata)
    return (avg, med, zmin, wheremin, zmax, wheremax)

def narrowdata(data, toprow, header):
    i = toprow.index(header)
    return [line[i] for line in data]


for file in os.listdir(datadir):
    print(file)
    infile = open(os.path.join(datadir, file), "r")
    header = infile.readline().strip().split("\t")
    lines = [destringify(line.strip().split("\t")) for line in infile]

    narrow = lambda x: narrowdata(lines, header, x)

    summary = {str:summarize(narrow(str)) for str in ['dn', 'ds', 'dnds', 'Pi']}


    outputfilename = file + ".summary"
    outfile = open(os.path.join("summary_stats", outputfilename), 'w')

    outfile.write("value\tmean\tmedian\tzmin\twheremin\tzmax\twheremax\n")
    for k,v in summary.items():
        outfile.write(k+'\t' + '\t'.join(map(str, v)) + '\n')
    outfile.close()
    infile.close()
    # dnvalues = [line[getindexof('dn')] for line in lines]
    # dsvalues = [line[getindexof('ds')] for line in lines]
    # dndsvalues = [line[getindexof('dnds')] for line in lines]
    # Pivalues = [line[getindexof('Pi')] for line in lines]



    
