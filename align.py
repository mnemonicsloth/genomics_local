

import os, sys

working_dir = 'ccout/Yersinia_pestis'



def process_dir(dir):
    os.mkdir(os.path.join(dir, 'align'))
    for file in os.listdir(os.path.join(dir, 'core')):
        corefile = os.path.join(dir, 'core', file)
        alignfile = os.path.join(dir, 'align', file + '.align')
        os.system('muscle -in ' + corefile + ' -out ' + alignfile)



def process_dir(dir):
    root = '/nas/longleaf/home/jdeszyck/pine'
    align_out = os.mkdir(os.path.join(root, 'align',  dir))
    for file in os.listdir(os.path.join( root, 'ccout', dir, 'core')):
        
        


