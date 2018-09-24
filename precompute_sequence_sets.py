"""
Precompute normalization and frequency vectors for public datasets, or any dataset using the FASTA format

@author: Chris Horning, Shuo Wang

"""
import os
import pathlib
import time

import numpy as np

import getGenCode
import skr_config
from seekr_launch_utils import compute_normalization_and_frequency
import pickle

CACHE_DIR = 'cache'
"""The file types that are persisted, and their respective name """
CACHE_FILE_TYPES = {'mean':'mean', 'std':'std', 'unnormalized_frequency':'unnormalized_frequency',
                    'names': 'names'}
NUMPY_EXTENSION = '.npy'
NAMES_EXTENSION = '.bin'

verbose = False

def initialize_cache():

    getGenCode.clear_cache()
    getGenCode.download_unzip(skr_config.GENCODE_HUMAN)
    getGenCode.download_unzip(skr_config.GENCODE_MOUSE)

    build_cache_files()

def get_precomputed_fasta_sets():
    fasta_sets = [skr_config.GENCODE_HUMAN, skr_config.GENCODE_MOUSE]
    return fasta_sets


def get_file_path_for(fasta_file, kmer_length, cache_file_type):
    dir_name = pathlib.PurePath(fasta_file).stem
    path_to_dir = os.path.join(CACHE_DIR, dir_name)

    if cache_file_type ==  CACHE_FILE_TYPES.get('names'):
        file_name = os.path.join(path_to_dir, cache_file_type + NAMES_EXTENSION)
    else:
        file_name = os.path.join(path_to_dir, cache_file_type + str(kmer_length) + NUMPY_EXTENSION)

    return file_name

def build_cache_files():
    if not os.path.exists(CACHE_DIR):
        os.mkdir(CACHE_DIR)

    fasta_sets = get_precomputed_fasta_sets()

    for fasta_set in fasta_sets:
        fasta_file = getGenCode.get_unzipped_file_name(fasta_set)
        dir_name = pathlib.PurePath(fasta_file).stem
        path_to_dir = os.path.join(CACHE_DIR, dir_name)
        if not os.path.exists(path_to_dir):
            os.mkdir(path_to_dir)

        names_written = False
        tsave = 0
        if verbose:
            print(dir_name + ' computing normalization took\t', end='')
        for kmer_length in range(1, skr_config.MAX_KMER_LENGTH_PRECOMPUTE + 1):
            fasta_path = os.path.join(CACHE_DIR, fasta_file)
            with open(fasta_path, mode='r') as infasta:
                t1 = time.perf_counter()
                (mean, std, unnormalized_frequency, names) = compute_normalization_and_frequency(infasta, kmer_length, return_normalized=False)
                t2 = time.perf_counter()
                if verbose:
                    print('k=' + str(kmer_length) + ',%.3fs;\t' % (t2 - t1), end='')

                t1 = time.perf_counter()
                np.save(get_file_path_for(fasta_file, kmer_length, CACHE_FILE_TYPES.get('mean')), mean)
                np.save(get_file_path_for(fasta_file, kmer_length, CACHE_FILE_TYPES.get('std')), std)
                np.save(get_file_path_for(fasta_file, kmer_length, CACHE_FILE_TYPES.get('unnormalized_frequency')), unnormalized_frequency)

                if not names_written:
                    with open(get_file_path_for(fasta_file, kmer_length, CACHE_FILE_TYPES.get('names')), 'wb') as names_file:
                        pickle.dump(names, names_file)
                    names_written = True

                t2 = time.perf_counter()
                tsave += (t2 - t1)
        if verbose:
            print('\nAggregate save time for ' + dir_name + ' was %.3fs' % tsave)

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        for arg in sys.argv[1:]:
            if arg == '-v':
                verbose = True
    t1 = time.perf_counter()
    initialize_cache()
    t2 = time.perf_counter()
    print('Initializing the cache took %.3f seconds' % (t2 - t1))

