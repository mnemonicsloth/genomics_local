import os
import os.path
import sys

in_dir = sys.argv[0]
out_dir = sys.argv[1]
ext = sys.argv[2]

def do_dir(dir):
    inpath = join(in_dir, dir)
    outpath = join(out_dir, dir)
    command =  f"python master.py -in {inpath} -out {outpath} -ext {ext}"
    os.system(command)


for f in os.listdir(in_dir):
    if os.path.isdir(f):
        do_dir(f)


