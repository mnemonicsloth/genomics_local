import os, sys


def count_lines(file):
   return sum(1 for line in open(file, 'r'))




small_cores = [(s, count_lines(os.path.join('ccout', s, 'families_core.txt'))) for s in os.listdir('ccout')
 if count_lines(os.path.join('ccout', s, 'families_core.txt')) < 500]
