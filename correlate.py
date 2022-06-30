import numpy as np
import scipy.stats as st
import os
import math

input_dir = 'pi_merge_output'
output_dir = 'corr'

def process_file(filename):
    infile = open(os.path.join(input_dir, filename), 'r')
    header = infile.readline().strip().split('\t')


    def get_index(column_name):
        return header.index(column_name)
        
    data = [ line.strip().split('\t') for line in infile ]

    def get_column(index):
        return [ line[index] for line in data ]

    def convert_column_to(f, column):
        return list(map(f, get_column(get_index(column))))

    avg_dnds = convert_column_to(float, 'avg_dnds')
    med_dnds = convert_column_to(float, 'med_dnds')
    dn = convert_column_to(float, 'dn')
    ds = convert_column_to(float, 'ds')
    Pi = convert_column_to(float, 'Pi')
    loc = convert_column_to(int, 'loc')

    def correlation(xs):
        pairs = [(x, l) for x, l in zip(xs, loc) if math.isnan(x)==False]
        newxs = [x for x,y in pairs]
        newloc = [y for x,y in pairs]
        return st.spearmanr(newxs, newloc)

    adndsr, adndsp = correlation(avg_dnds)
    mdndsr, mdndsp = correlation(med_dnds)
    dnr, dnp = correlation(dn)
    dsr, dsp = correlation(ds)
    pir, pip = correlation(Pi)

    infile.close()
    return adndsr, adndsp, mdndsr, mdndsp, dnr, dnp, dsr, dsp, pir, pip



def go():
    outfile = open(os.path.join(output_dir, 'correlations.txt'), 'w')
    outfile.write('avg_dndsr\tavg_dndsp\tmed_dndsr\tmed_dndsp\tdnr\tdnp\tdsr\tdsp\tpir\tpip\n')

    for file in os.listdir(input_dir):
        print(file)
        adndsr, adndsp, mdndsr, mdndsp, dnr, dnp, dsr, dsp, pir, pip = process_file(file)
        outfile.write(f'{file}\t{adndsr}\t{adndsp}\t{mdndsr}\t{mdndsp}\t{dnr}\t{dnp}\t{dsr}\t{dsp}\t{pir}\t{pip}\n')

    outfile.close()

    



        
