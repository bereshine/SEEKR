# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 23:48:48 2017

@author: Chris
"""

import os
from io import BytesIO
from io import StringIO

import numpy as np

import getGenCode
import kmer_counts
import skr_config
from SeekrServerError import SeekrServerError
from precompute_sequence_sets import CACHE_FILE_TYPES
from precompute_sequence_sets import get_file_path_for
from precompute_sequence_sets import get_precomputed_fasta_sets
from seekr_launch_utils import compute_normalization_and_frequency
from seekr_launch_utils import get_names_from_counter
from pearson import pearson
import pandas
import pickle
from session_helper import get_file_for_directory_id
from io import TextIOWrapper

def run_seekr_algorithm(parameters):
    """
    Launch the SEEKR algorithm using the parameters from the Web Service and return a zip file of the results

    The logic in this method tries to avoid recalculating k-mers and normalization for performance reasons.  The logic
    otherwise would be simpler -
    1) Calculate the normalization (mean and std dev)
    2) Calculate the frequencies of the user set and then apply the normalization from step 1
    3) Calculate the frequencies of the comparison set if it exists and apply the normalization from step 1
    4) Calculate the Pearson's R correlations between sequences in the user set and the comparison set.
        If no comparison set exists, calculate the correlations between the user set sequences

    In any of these steps, if we already have a precomputed value, we will load that instead of performing the computation.

    Notes
    -----
    numpy's corrcoef is an efficient way to calculate Pearson's correlations, but since its implementation computes a
    covariance matrix, the output is always a square matrix.  So if we had 10 sequences in a user set and compare
    against 10,000 sequences in a comparision set, numpy.corrcoef will calculate a matrix that is 10,010x10,010.
    The pearson function called supports non-square matrices and is thus used for comparing against the comparision set.
    e.g. it's matrix would be 10x10,000.

    """
    outfile = 'test1.csv'
    mean_std_loaded = False
    names = None
    comparison_names = None
    normalization_path = get_precomputed_normalization_path(parameters)
    if normalization_path is not None:
        mean = np.load(normalization_path[0])
        std = np.load(normalization_path[1])
        mean_std_loaded = True

    normal_set = parameters['normal_set']
    if normal_set is None:
        raise SeekrServerError('No normalization set Provided')
    comparison_set = None
    if 'comparison_set' in parameters:
        comparison_set = parameters['comparison_set']
    if 'comparison_set_files' in parameters:
        if normal_set == skr_config.SETTING_USER_SET:
            (mean, std, counts, names) = compute_normalization_and_frequency(
                infasta=TextIOWrapper(parameters['user_set_files']), kmer_length=parameters['kmer_length'], outfile=outfile)
            counter = kmer_counts.BasicCounter(infasta=TextIOWrapper(parameters['comparison_set_files']), outfile=None,
                                               k=parameters['kmer_length'],
                                               label=True, silent=True, binary=False, mean=mean, std=std)
            comparison_counts = counter.make_count_file()
            comparison_names = get_names_from_counter(counter)
        elif normal_set == skr_config.SETTING_COMPARISION_SET:
            (mean, std, comparison_counts, comparison_names) = compute_normalization_and_frequency(
                infasta=TextIOWrapper(parameters['comparison_set_files']), kmer_length=parameters['kmer_length'])
            counter = kmer_counts.BasicCounter(infasta=TextIOWrapper(parameters['user_set_files']), outfile=outfile,
                                               k=parameters['kmer_length'],
                                               label=True, silent=True, binary=False, mean=mean, std=std)
            counts = counter.make_count_file()
            names = get_names_from_counter(counter)

        elif mean_std_loaded:
            counter = kmer_counts.BasicCounter(infasta=TextIOWrapper(parameters['user_set_files']), outfile=outfile,
                                               k=parameters['kmer_length'],
                                               label=True, silent=True, binary=False, mean=mean, std=std)
            counts = counter.make_count_file()

            comparision_counter = kmer_counts.BasicCounter(infasta=TextIOWrapper(parameters['comparison_set_files']), outfile=None,
                                               k=parameters['kmer_length'],
                                               label=True, silent=True, binary=False, mean=mean, std=std)
            comparison_counts = comparision_counter.make_count_file()
            names = get_names_from_counter(counter)
            comparison_names = get_names_from_counter(comparision_counter)

        else:
            raise SeekrServerError('Normalization for Comparision Set File is not valid')

        similarity = pearson(counts, comparison_counts)
    elif comparison_set is not None and len(comparison_set) > 0 and comparison_set != 'user_set':
        unnormalized_frequency_path, names_path = get_precomputed_frequency_path(comparison_set, parameters['kmer_length'])
        assert unnormalized_frequency_path is not None and names_path is not None

        if normal_set == skr_config.SETTING_USER_SET:
            (mean, std, counts, names) = compute_normalization_and_frequency(
                infasta=TextIOWrapper(parameters['user_set_files']), kmer_length=parameters['kmer_length'], outfile=outfile)
            counter = kmer_counts.BasicCounter(infasta=TextIOWrapper(parameters['comparison_set_files']), outfile=None,
                                               k=parameters['kmer_length'],
                                               label=True, silent=True, binary=False, mean=mean, std=std)
            comparison_counts = _unnormalized_frequency_to_normalized(unnormalized_frequency_path, mean, std)
            comparison_names = load_names_from_path(names_path)
        elif normal_set == skr_config.SETTING_COMPARISION_SET:
            raise SeekrServerError('')

        elif mean_std_loaded:
            counter = kmer_counts.BasicCounter(infasta=TextIOWrapper(parameters['user_set_files']), outfile=outfile,
                                               k=parameters['kmer_length'],
                                               label=True, silent=True, binary=False, mean=mean, std=std)
            counts = counter.make_count_file()

            comparison_counts = _unnormalized_frequency_to_normalized(unnormalized_frequency_path, mean, std)
            names = get_names_from_counter(counter)
            comparison_names = load_names_from_path(names_path)

        else:
            raise SeekrServerError('No normalization set Provided')

        similarity = pearson(counts, comparison_counts)

    else:
        if mean_std_loaded:
            counter = kmer_counts.BasicCounter(infasta=TextIOWrapper(parameters['user_set_files']), outfile=outfile,
                                               k=parameters['kmer_length'],
                                               label=True, silent=True, binary=False, mean=mean, std=std)
            counts = counter.make_count_file()
        elif normal_set == skr_config.SETTING_USER_SET:
            counter = kmer_counts.BasicCounter(infasta=TextIOWrapper(parameters['user_set_files']), outfile=outfile, k=parameters['kmer_length'],
                                           label=True, silent=True, binary=False)
            counts = counter.make_count_file()
        else:
            raise SeekrServerError('Normalization type is not valid')

        names = get_names_from_counter(counter)
        similarity = np.corrcoef(counts)

    #TODO refactor - original code saved to csv on disk - move this to a separate operation
    with open(outfile) as csvFile:
        counts_text = csvFile.read()

    bytes_io = BytesIO()
    np.save(bytes_io, similarity)
    bytes_io.seek(0)
    pearsons_file_in_memory = bytes_io.read()

    return counts_text, pearsons_file_in_memory

def _run_seekr_algorithm(parameters):
    """
    Launch the SEEKR algorithm using the parameters from the Web Service and return a zip file of the results

    The logic in this method tries to avoid recalculating k-mers and normalization for performance reasons.  The logic
    otherwise would be simpler -
    1) Calculate the normalization (mean and std dev)
    2) Calculate the frequencies of the user set and then apply the normalization from step 1
    3) Calculate the frequencies of the comparison set if it exists and apply the normalization from step 1
    4) Calculate the Pearson's R correlations between sequences in the user set and the comparison set.
        If no comparison set exists, calculate the correlations between the user set sequences

    In any of these steps, if we already have a precomputed value, we will load that instead of performing the computation.

    Notes
    -----
    numpy's corrcoef is an efficient way to calculate Pearson's correlations, but since its implementation computes a
    covariance matrix, the output is always a square matrix.  So if we had 10 sequences in a user set and compare
    against 10,000 sequences in a comparision set, numpy.corrcoef will calculate a matrix that is 10,010x10,010.
    The pearson function called supports non-square matrices and is thus used for comparing against the comparision set.
    e.g. it's matrix would be 10x10,000.

    """
    outfile = None
    mean_std_loaded = False
    names = None
    comparison_counts = None
    comparison_names = None
    normalization_path = get_precomputed_normalization_path(parameters)
    if normalization_path is not None:
        mean = np.load(normalization_path[0])
        std = np.load(normalization_path[1])
        mean_std_loaded = True

    normal_set = parameters['normal_set']
    if normal_set is None:
        raise SeekrServerError('No normalization set Provided')
    comparison_set = None
    if 'comparison_set' in parameters:
        comparison_set = parameters['comparison_set']
    if 'comparison_set_files' in parameters:
        if normal_set == skr_config.SETTING_USER_SET:
            (mean, std, counts, names) = compute_normalization_and_frequency(
                infasta=_load_user_set_file(parameters), kmer_length=parameters['kmer_length'], outfile=outfile)
            counter = kmer_counts.BasicCounter(infasta=_load_comparison_set_file(parameters), outfile=None,
                                               k=parameters['kmer_length'],
                                               label=True, silent=True, binary=False, mean=mean, std=std)
            comparison_counts = counter.make_count_file()
            comparison_names = get_names_from_counter(counter)
        elif normal_set == skr_config.SETTING_COMPARISION_SET:
            (mean, std, comparison_counts, comparison_names) = compute_normalization_and_frequency(
                infasta=_load_comparison_set_file(parameters), kmer_length=parameters['kmer_length'])
            counter = kmer_counts.BasicCounter(infasta=_load_user_set_file(parameters), outfile=outfile,
                                               k=parameters['kmer_length'],
                                               label=True, silent=True, binary=False, mean=mean, std=std)
            counts = counter.make_count_file()
            names = get_names_from_counter(counter)

        elif mean_std_loaded:
            counter = kmer_counts.BasicCounter(infasta=_load_user_set_file(parameters), outfile=outfile,
                                               k=parameters['kmer_length'],
                                               label=True, silent=True, binary=False, mean=mean, std=std)
            counts = counter.make_count_file()

            comparision_counter = kmer_counts.BasicCounter(infasta=_load_comparison_set_file(parameters), outfile=None,
                                               k=parameters['kmer_length'],
                                               label=True, silent=True, binary=False, mean=mean, std=std)
            comparison_counts = comparision_counter.make_count_file()
            names = get_names_from_counter(counter)
            comparison_names = get_names_from_counter(comparision_counter)

        else:
            raise SeekrServerError('Normalization for Comparision Set File is not valid')

    elif comparison_set is not None and len(comparison_set) > 0 and comparison_set != 'user_set':

        unnormalized_frequency_path, names_path = get_precomputed_frequency_path(comparison_set, parameters['kmer_length'])
        assert unnormalized_frequency_path is not None and names_path is not None

        if normal_set == skr_config.SETTING_USER_SET:
            (mean, std, counts, names) = compute_normalization_and_frequency(
                infasta=_load_user_set_file(parameters), kmer_length=parameters['kmer_length'], outfile=outfile)
            counter = kmer_counts.BasicCounter(infasta=_load_comparison_set_file(parameters), outfile=None,
                                               k=parameters['kmer_length'],
                                               label=True, silent=True, binary=False, mean=mean, std=std)
            comparison_counts = _unnormalized_frequency_to_normalized(unnormalized_frequency_path, mean, std)
            comparison_names = load_names_from_path(names_path)
        elif normal_set == skr_config.SETTING_COMPARISION_SET:
            raise SeekrServerError('')

        elif mean_std_loaded:
            counter = kmer_counts.BasicCounter(infasta=_load_user_set_file(parameters), outfile=outfile,
                                               k=parameters['kmer_length'],
                                               label=True, silent=True, binary=False, mean=mean, std=std)
            counts = counter.make_count_file()

            comparison_counts = _unnormalized_frequency_to_normalized(unnormalized_frequency_path, mean, std)
            names = get_names_from_counter(counter)
            comparison_names = load_names_from_path(names_path)

        else:
            raise SeekrServerError('No normalization set Provided')

    else:
        if mean_std_loaded:
            counter = kmer_counts.BasicCounter(infasta=_load_user_set_file(parameters), outfile=outfile,
                                               k=parameters['kmer_length'],
                                               label=True, silent=True, binary=False, mean=mean, std=std)
            counts = counter.make_count_file()
        elif normal_set == skr_config.SETTING_USER_SET:
            counter = kmer_counts.BasicCounter(infasta=_load_user_set_file(parameters), outfile=outfile, k=parameters['kmer_length'],
                                           label=True, silent=True, binary=False)
            counts = counter.make_count_file()
        else:
            raise SeekrServerError('Normalization type is not valid')

        names = get_names_from_counter(counter)

    return counts, names, comparison_counts, comparison_names, counter

def get_precomputed_normalization_path(parameters):
    normal_set = parameters['normal_set']

    if normal_set is None or len(normal_set) <= 0:
        return None

    fasta_sets =  get_precomputed_fasta_sets()
    for fasta_set in fasta_sets:
        if normal_set == fasta_set.server_name:
            fasta_file = getGenCode.get_unzipped_file_name(fasta_set)
            mean_path = get_file_path_for(fasta_file, parameters['kmer_length'], CACHE_FILE_TYPES.get('mean'))
            std_path = get_file_path_for(fasta_file, parameters['kmer_length'], CACHE_FILE_TYPES.get('std'))
            if os.path.exists(mean_path) and os.path.exists(std_path):
                return (mean_path, std_path)
            else:
                raise SeekrServerError('Fasta file <' + fasta_file + '> not found for kmer_length=' +  str(parameters['kmer_length']) )
                return None

    return None

def get_precomputed_frequency_path(comparison_set, kmer_length):
    if comparison_set is None or len(comparison_set) <= 0:
        return None

    fasta_sets =  get_precomputed_fasta_sets()
    for fasta_set in fasta_sets:
        if comparison_set == fasta_set.server_name:
            fasta_file = getGenCode.get_unzipped_file_name(fasta_set)
            unnormalized_frequency_path = get_file_path_for(fasta_file, kmer_length, CACHE_FILE_TYPES.get('unnormalized_frequency'))
            names_path = get_file_path_for(fasta_file, kmer_length, CACHE_FILE_TYPES.get('names'))
            if os.path.exists(unnormalized_frequency_path) and os.path.exists(names_path):
                return unnormalized_frequency_path, names_path

    return None

def _unnormalized_frequency_to_normalized(unnormalized_frequency_path, mean, std):
    normalized_frequency = np.load(unnormalized_frequency_path)
    normalized_frequency -= mean
    normalized_frequency /= std
    return normalized_frequency

def get_data_for_pearsons(counts, comparison_counts, col1_names, col2_names):
    similarity = pearson(counts, comparison_counts)
    df = pandas.DataFrame(data=similarity, index=col1_names, columns=col2_names)
    return df

def load_names_from_path(names_path):
    with open(names_path, 'rb') as names_file:
        names = pickle.load(names_file)
    return names

def _load_user_set_file(parameters):
    file = get_file_for_directory_id(parameters['directory_id'], parameters['user_set_files'], extension='fasta')
    return StringIO(file)

def _load_comparison_set_file(parameters):
    file = get_file_for_directory_id(parameters['directory_id'], parameters['comparison_set_files'], extension='fasta')
    return StringIO(file)


def fixup_counts(counts, counter):
    """
    Set all values where std=0 or NaN to zero.
    """
    warnings = list()
    assert counter.std is not None

    zero_indexes = np.where(counter.std == 0)
    if zero_indexes[0].size > 0:
        counts[:, zero_indexes] = 0
        warnings.append('Warning: Normalization contains k-mers with 0 standard deviation')

    nan_indexes = np.where(np.isnan(counter.std))
    if nan_indexes[0].size > 0:
        counts[:, nan_indexes] = 0
        warnings.append('Warning: Normalization contains k-mers with NaN standard deviation')

    return warnings

def get_kmers_csv(counts, names, kmers):
    df = pandas.DataFrame(data=counts, index=names, columns=kmers)
    return df.to_csv()

def get_pearsons_csv(names, pearsons, comparison_names):
    df = pandas.DataFrame(data=pearsons, index=names, columns=comparison_names)
    return df.to_csv()