#! /usr/bin/python
import os
import sys
from os import listdir
from os.path import isfile, join, basename
import re
import shutil


fadir = "results.fa"
fa_stems = [ re.sub('_genomic\.fa$', '', f) for f in listdir(fadir) if isfile(join(fadir, f))]


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


def make_fa_filename(fa_stem):
    return fa_stem + '_genomic.fa'

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




fa_dir = "/home/jd/src/file-of-genomes/results.fa"
output_dir = "/home/jd/src/file-of-genomes/species/"
def make_speciated_files(dict):
    for species, entries in abridged_dict.items():
        species_dir = os.path.join(output_dir, species)
        os.mkdir(species_dir)
        for stem in entries:
            fname = make_fa_filename(stem)
            shutil.copyfile(os.path.join(fa_dir, fname), os.path.join(species_dir, fname))


abridged_dict = abridge(make_species_dictionary(fa_stems), 5)

def main():
    make_speciated_files(abridged_dict)


        















