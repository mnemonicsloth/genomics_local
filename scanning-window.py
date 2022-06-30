import numpy as np

datadir = "pi_merge_output"
windowsize = 50000

def getIndexOf(str):
    return header.index(str)


def destringify(line):
    tofloats = [1,2,3,4,5,6,7,8,9,11]
    for i in tofloats:
        line[i] = float(line[i])
    line[locidx] = int(line[locidx])
    return line

def inrange(n, bounds):
    low, high = bounds
    if low <= n <= high:
        return True
    else:
        return False

def buildranges(maximum):
    step = int(windowsize / 2)
    lowers = range(0, maximum, step)
    uppers = range(windowsize, maximum + windowsize, step)
    return list(zip(lowers, uppers))

def findranges(n, ranges):
    return [bounds for bounds in ranges if inrange(n, bounds)]


def tallyrange(lines, property):
    index = getIndexOf(property)
    data = [line[index] for line in lines]
    if len(data) == 0:
        return 0.0
    else: 
        return sum(data)/ len(data)





for file in os.listdir(datadir):
    infile = open(os.path.join(datadir, file), "r")
    header = infile.readline().strip().split("\t")
    locidx = getIndexOf('loc')



    lines = [destringify(line.strip().split("\t")) for line in infile]
    inorderlines = sorted(lines, key=lambda line: line[10])

    lastloc = inorderlines[-1][locidx]

    ranges = buildranges(lastloc)
    bins = {r:[] for r in ranges}
    for line in lines:
        n = line[locidx]
        for r in findranges(n, ranges):
            bins[r].append(line)

    results = {}
    for rng, lines in bins.items():
        avgdnds = tallyrange(lines, "avg_dnds")
        dn  = tallyrange(lines, "dn")
        ds = tallyrange(lines, "ds")
        Pi = tallyrange(lines, "Pi")
        results[rng] = (avgdnds, dn, ds, Pi)

    ranges = sorted(list(results.keys()))


    outfile = open(os.path.join("sliding_window", file+".window" ), "w")
    outfile.write(f"low-high\tdnds\tdn\tds\tPi\n")
    for low,high in ranges:
        dnds, dn, ds, Pi = results[(low,high)]
        outfile.write(f"{low}-{high}\t{dnds}\t{dn}\t{ds}\t{Pi}\n")
    outfile.close()

