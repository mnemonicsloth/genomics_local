import os, sys

basedir = '/home/jd/src/file-of-genomes'
genomedir = os.path.join(basedir, 'species.prot')
ccoutdir = os.path.join(basedir, 'ccout')
outfile = os.path.join(basedir, 'core-genomes.txt')
of = open(outfile, 'w')



for f in os.listdir(genomedir):
    species = f
    num_genomes = len(os.listdir(os.path.join(genomedir, f)))
    core_dir = os.path.join(ccoutdir, f)
    corefile = os.path.join(core_dir, 'families_core.txt')
    if not os.path.exists(corefile):
        count = 0
    else:
        count = len(open(corefile).readlines())
    of.write(f"{species:40s}{num_genomes:10d}{count:10d}\n")

    
of.close()

    
