import os
import math
import scipy.stats as ss

input_dir = 'sliding_window'
output_dir = 'zscores.reduced'
num_lines = 10

def get_extremal_values(xs, n, direction):
    topN = []
    for i in range(len(xs)):
        topN.append( (xs[i], i) )
        if direction == 'highest':
            topN.sort(reverse=True)
        else:
            topN.sort()
        if len(topN) > n:
            topN = topN[:-1]
    return topN


def process_file(filename):
    infile = open(os.path.join(input_dir, filename), 'r')
    header = infile.readline().strip().split('\t')

    def get_index(colname):
        return header.index(colname)

    data = [line.strip().split('\t') for line in infile]
    for line in data:
        for i in range(1, len(line)):
            line[i] = float(line[i])

    data = [line for line in data if not any(map(math.isnan, line[1:]))]
    if len(data) == 0:
        return None

    def docolumn(colname):
        return [line[get_index(colname)] for line in data]

    cols = list(map(docolumn, header))
    datadict = dict(zip(header, cols))

    for key in header[1:]:
        datadict[key+'-z'] = list(ss.zscore(datadict[key], nan_policy="omit"))

    for key in header[1:]:
        pvals = [ss.norm.sf(x)*2 for x in datadict[key+'-z']]
        datadict[key+'-p'] = pvals


    datadict['low-high'] = [int(x.split('-')[0]) for x in datadict['low-high']]


    
    long_header = ['low-high', 'dnds', 'dn', 'ds', 'Pi', 'dnds-z', 'dn-z', 'ds-z', 'Pi-z', 'dnds-p', 'dn-p', 'ds-p', 'Pi-p']
        
    def get_row(i, d):
        return [val[i] for key,val in d.items()]

    def get_column_extremal_values(colname, df):
        maxes = get_extremal_values(df[colname], num_lines, 'highest')
        mins =  get_extremal_values(df[colname], num_lines, 'lowest')

        maxes = [get_row(i, df) for (x,i) in maxes]
        mins =  [get_row(i, df) for (x,i) in mins]

        return maxes, mins

    def write_column_output(colname):
        highs, lows = get_column_extremal_values(colname, datadict)
        header_out = '\t'.join(long_header) + '\n'

        output_stem = os.path.join(output_dir, filename + '.' + colname)
        high_filename = output_stem + '.' + 'high'
        low_filename =  output_stem + '.' + 'low'

        outfile = open(high_filename, 'w')
        print(high_filename)
        outfile.write(header_out)
        for line in highs:
            line = map(str, line)
            outfile.write('\t'.join(line) + '\n')
        outfile.close()

        outfile = open(low_filename, 'w')
        print(low_filename)
        outfile.write(header_out)
        for line in lows:
            line = map(str, line)
            outfile.write('\t'.join(line) + '\n')
        outfile.close()
        return header_out
        
    for colname in header[1:]:
        write_column_output(colname)
    

    return 1


def main():
    for file in os.listdir(input_dir):
        process_file(file)
            


    
