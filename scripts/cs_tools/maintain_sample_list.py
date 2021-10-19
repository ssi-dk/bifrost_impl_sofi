#!/usr/bin/env python3

import argparse
from pathlib import Path
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
    _sample_list_path: Path

    def __init__(self, sample_list_path: Path):
        self._sample_list_path = sample_list_path
        if self._sample_list_path.exists():
            print(f"Reading sample list in {str(self._sample_list_path)}.")
            with open(self._sample_list_path, 'r') as sample_list:
                for (sample_name, file1, file2) in read_sample_list(sample_list):
                    if sample_name == 'sample':  # ignore header
                        continue
                    self.add_sample(sample_name, file1, file2)
        else:
            print(f"{str(self._sample_list_path)} did not exist, but will be created.")
    
    def add_sample(self, sample_name: str, file1: str, file2: str):
        # if sample_name in self._samples:
        #    print(f"{sample_name} is already in sample list.")
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
            if len(line) > 1:  # If not in last line (which should only contain a '\n')
                sample_name, file1, file2 = line.rstrip().split('\t')
                output.append((sample_name, file1, file2))
        except StopIteration:
            break
    return output

def find_new_samples(fastq_dir: Path, metadata_file: Path, species: str):
    """
    Find samples of provided species in folder fastq_dir.
    This is done by, looking up species in run_metadata.tsv and filter by 'species' in second column.
    Return a list of (sample_name, file1, file2) where file1, file2 have full paths.
    """
    output = list()
    r1_files_gen = fastq_dir.glob("*_S*_R1_001.fastq*")
    with open(metadata_file, 'r') as run_metadata_tsv:
        species_lookup_dict = dict()
        for line in run_metadata_tsv:
            sample_species = line.split('\t')[0:2]
            try:
                species_lookup_dict[sample_species[0]] = sample_species[1]
            except IndexError:
                print(f"Error when looking for samples in fastq_dir {fastq_dir}")
                print("sample_species list:")
                print(sample_species)
    
    while True:
        try:
            file1_path = Path(next(r1_files_gen))
            sample_name = file1_path.stem[:-22]
            try:
                found_species: str = species_lookup_dict[sample_name.replace("'", "")]
            except KeyError:
                print(f"For the sample {file1_path} we could not look up a species in the metadata file.")
                found_species = "Unknown"
            if not found_species == species:
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

def maintain(species: str, sample_list: Path, fastq_dir: Path) -> bool:
    metadata_file = Path(fastq_dir.joinpath('run_metadata.tsv'))
    if not metadata_file.exists():
        print(f"Metadata file {metadata_file} does not exist - ignoring folder {fastq_dir}")
        return False
    container = SampleContainer(sample_list)
    fastq_dir = Path(fastq_dir)
    species = str(species).replace('_', ' ')
    new_samples = find_new_samples(fastq_dir, metadata_file, species)
    if len(new_samples) == 0:
        print(f"No new fastq files found in folder {fastq_dir} - {sample_list} not modified.")
    else:
        print("NEW SAMPLES:")
        print_samples(new_samples)
        for new_sample in new_samples:
            container.add_sample(*new_sample)
        container.save()
    return True

# Running maintain_sample_list.py as a standalone script is currently not maintained.
""" if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Create and maintain a global sample list for chewieSnake.")
    parser.add_argument('sample_list', help="Path and filename for sample list. If file does not exist it will be created.")
    parser.add_argument('-d', '--fastq_dir', help="Path to existing directory containing fastq files. Default: current directory.")
    parser.add_argument('-s', '--species', help="Optional species name to filter on. Spaces between name components must be \
        replaced with underscores; f. ex. Salmonella_enterica.\
        If provided, a lookup will be made in a file named 'run_metadata.tsv' inside fastq_dir, and only those \
        samples where --species - with underscores replaced by spaces - match the second column in this file will be added to the \
        sample list.")
    args = parser.parse_args()
    sample_list = Path(args.sample_list)
    fastq_dir = Path(args.fastq_dir)
    maintain(args.species, sample_list, fastq_dir)
 """
