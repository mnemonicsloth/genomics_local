import os
import math
import scipy.stats as ss
import pandas as pd

loc_source = "pi_merge_output"
zscore_source = "zscores.reduced"
output_dir = "zscores.reduced.sorted"



def get_short_filename(filename):
    fn = filename.split('.')
    return fn[0] + '.txt'

test_fn = 'Escherichia_coli.txt.window.dnds.high'

def get_max_gene_loc(filename):
    shortfn = get_short_filename(filename)
    shortfn = os.path.join(loc_source, shortfn)
    df = pd.read_csv(shortfn, sep='\t', header=0)
    return df['loc'].max()

#    pd.read_csv('', sep='\t', header=0)
