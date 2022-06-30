import os, sys


sfa = 'species.fa'
sprot = 'species.prot'

for dir in os.listdir('species.prot'):
    newpath = os.path.join('ccout', dir, 'species.prot')
    # os.mkdir(newpath)
    for file in os.listdir(os.path.join('species.prot', dir)):
        os.rename(os.path.join('species.prot', dir, file),
                  os.path.join(newpath, file))

