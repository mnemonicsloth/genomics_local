import os
import sys


def pick_reference_genome(directory):
    return os.path.join(directory, os.listdir(directory)[0])


def get_gene_ids(file):
    base = os.path.basename(file)
    return  [ base + '&' + line.strip(">- \n") \
              for line in open(file, 'r') \
              if re.search("^>", line)]


def process_families_line(line, dict):
    genes = line.strip('\n').split('\t')
    for g in genes[1:]:
        if g in dict:
            raise Exception("gene assigned to multiple IDs")
        # print(f"Adding {g} to dict")
        dict[g] = genes[0]
    return dict


def process_families_file(file):
    dict = {}
    for line in open(file, 'r'):
        dict = process_families_line(line, dict)
    return dict


def do_species(dir):
    genes = get_gene_ids(pick_reference_genome(os.path.join(dir, 'species.prot')))
    dict = process_families_file(os.path.join(dir, 'families_core.txt'))

    of = open(os.path.join(dir, 'ref_core.txt'), 'w')

    i = 0
    for gene in genes:
        if gene in dict:
            of.write("{:<10} {:<10} {}\n".format(dict[gene], i, gene))
        i += 1


def main():
    for species in os.listdir('ccout'):
        print(species + '\n')
        do_species(os.path.join('ccout', species))
