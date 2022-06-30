#! /usr/bin/python
import os
import sys
from os import listdir
from os.path import isfile, join, basename
import re
import shutil


prdir = "results.prot"
prot_stems = [ re.sub('_genomic\.prot$', '', f) for f in listdir(prdir) if isfile(join(prdir, f))]


record_file = "completes.csv"
records = [line for line in open(record_file, "r")]

def find_record(stem):
    pattern = re.compile(stem)
    matches = []
    for line in records:
        if pattern.search(line):
            matches.append(line)
    if len(matches) == 0:
        raise Exception(stem + ": no matches found!")
    elif len(matches) > 1:
        raise Exception(stem + ": " + "found " + len(matches) + " matches, not 1")
    else:
        return matches[0]

def read_species(record):
    words = re.match('".*?"', record).group().strip('"').split(" ")
    return words[0] + "_" + words[1]


def make_prot_filename(prot_stem):
    return prot_stem + '_genomic.prot'

def make_species_dictionary(stems):
    dict = {}
    for s in stems:
        species = read_species(find_record(s))
        if species in dict.keys():
            dict[species].append(s)
        else:
            dict[species] = [s]
    return dict


def abridge(dict, min_population):
    return {k:v for k,v in dict.items() if len(v) >= min_population}


# I miss reduce
def count_entries(dict):
    count = 0
    for file_list in dict.values():
        count += len(file_list)
    return count




prot_dir = "/home/jd/src/file-of-genomes/results.prot"
output_dir = "/home/jd/src/file-of-genomes/species.prot/"
def make_speciated_files(dict):
    for species, entries in abridged_dict.items():
        species_dir = os.path.join(output_dir, species)
        os.mkdir(species_dir)
        for stem in entries:
            pname = make_prot_filename(stem)
            shutil.copyfile(os.path.join(prot_dir, pname), os.path.join(species_dir, pname))

 
abridged_dict = abridge(make_species_dictionary(prot_stems), 5)

def main():
    make_speciated_files(abridged_dict)


        



















