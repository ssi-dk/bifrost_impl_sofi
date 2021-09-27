#!/usr/bin/env python3

import argparse
import pathlib
import os
import sys
from typing import Iterable

"""
Create and maintain a sample list for use with chewieSnake.

Currently only fastq files are supported, and fastq file names must follow the
Illumina standard:
    
    {samplename}_S*_R{1,2}_001.fastq*

This script only requires Python 3 with no additional packages.

Author: Finn Gruwier Larsen, figl@ssi.dk
"""


class SampleContainer:
    _samples: dict = {}
    _sample_list_path: pathlib.Path

    def __init__(self, sample_list_path: pathlib.Path):
        self._sample_list_path = sample_list_path
        if self._sample_list_path.exists():
            print(f"{str(self._sample_list_path)} exists - reading file.")
            with open(self._sample_list_path, 'r') as sample_list:
                for (sample_name, file1, file2) in read_sample_list(sample_list):
                    if sample_name == 'sample':  # ignore header
                        continue
                    self.add_sample(sample_name, file1, file2)
        else:
            print(f"{str(self._sample_list_path)} does not exist - new file.")
    
    def add_sample(self, sample_name: str, file1: str, file2: str):
        if sample_name in self._samples:
            print(f"Omitting {sample_name} as it is already in sample list.")
        self._samples[sample_name] = (file1, file2)
    
    def list_samples(self):
        items = ((k, v[0], v[1]) for k, v in self._samples.items())
        return items

    def save(self):
        with open(self._sample_list_path, 'w') as sample_list:
            sample_list.write('sample\tfq1\tfq2\n')
            for k, v in self._samples.items():
                line = '\t'.join((k, v[0], v[1])) + '\n'
                sample_list.write(line)
            sample_list.write('\n')
            

def read_sample_list(sample_list):
    output = list()
    while True:
        try:
            line: str = next(sample_list)
            sample_name, file1, file2 = line.rstrip().split('\t')
            output.append((sample_name, file1, file2))
        except StopIteration:
            break
    return output

def find_new_samples(fastq_dir: pathlib.Path, species: str = None):
    """
    Find samples in folder fastq_dir.
    If spec is provided, lookup species in run_metadata.tsv and filter by spec in second column.
    Return a list of (sample_name, file1, file2) where file1, file2 have full paths.
    """
    r1_files_gen = fastq_dir.glob("*_S*_R1_001.fastq*")
    if species:
        with open(fastq_dir.joinpath('run_metadata.tsv'), 'r') as run_metadata_tsv:
            species_lookup_dict = dict()
            for line in run_metadata_tsv:
                sample_species = line.split('\t')[0:2]
                species_lookup_dict[sample_species[0]] = sample_species[1]
    output = list()
    while True:
        try:
            file1_path = pathlib.Path(next(r1_files_gen))
            sample_name = file1_path.stem[:-22]
            if species:
                found_species: str = species_lookup_dict[sample_name.replace("'", "")]
                if not found_species == species:
                    print(f"Ignoring {sample_name} - species '{found_species}' does not match '{species}'.")
                    continue
            file1_filename = file1_path.parts[-1]
            listified_filename = list(file1_filename)
            listified_filename[-14] = '2'
            file2_filename = ''.join(listified_filename)
            file2_path = file1_path.with_name(file2_filename)
            assert(file2_path.exists())
            output.append((sample_name, str(file1_path), str(file2_path)))
        except StopIteration:
            break
    return output

def print_samples(sample_iter: Iterable):
    print("(SampleID , File 1)")
    for sample in sample_iter:
        print(sample[0], sample[1])

def main():
    parser = argparse.ArgumentParser(description="Create and maintain a global sample list for chewieSnake.")
    parser.add_argument('sample_list', help="Path and filename for sample list. If file does not exist it will be created.")
    parser.add_argument('-d', '--fastq_dir', help="Path to existing directory containing fastq files. Default: current directory.")
    parser.add_argument('-s', '--species', help="Optional species name to filter on. Spaces between name components must be \
        replaced with underscores; f. ex. Salmonella_enterica.\
        If provided, a lookup will be made in a file named 'run_metadata.tsv' inside fastq_dir, and only those \
        samples where --species - with underscores replaced by spaces - match the second column in this file will be added to the \
        sample list.")
    args = parser.parse_args()
    sample_list_path = pathlib.Path(args.sample_list)
    print(f"Sample list path: {sample_list_path}")
    container = SampleContainer(sample_list_path)
    print("OLD SAMPLES:")
    print_samples(container.list_samples())
    fastq_dir = pathlib.Path(args.fastq_dir or os.getcwd())
    print(f"Folder to add fastq files from: {fastq_dir}")
    if args.species:
        species = str(args.species).replace('_', ' ')
        print(f"Filter on species name: {species}.")
    else:
        print("Do not filter on species name.")
    new_samples = find_new_samples(fastq_dir, species)
    if len(new_samples) == 0:
        print(f"No fastq files found in folder {fastq_dir} - {sample_list_path} not modified.")
    else:
        print("NEW SAMPLES:")
        print_samples(new_samples)
        for new_sample in new_samples:
            container.add_sample(*new_sample)
        container.save()

if __name__ == '__main__':
    main()
