
"""
@author: Chris Horning
"""
import kmer_counts


def compute_normalization(infasta, kmer_length):
    counter = kmer_counts.BasicCounter(infasta, outfile=None,
                                       k=kmer_length, label=True, silent=True, binary=True,
                                       mean=True, std=True)
    counter.make_count_file()
    return (counter.mean, counter.std)


def compute_normalization_and_frequency(infasta, kmer_length, return_normalized=True, outfile=None):
    counter = kmer_counts.BasicCounter(infasta, outfile=outfile,
                                       k=kmer_length, label=True, silent=True, binary=False,
                                       mean=True, std=True)
    counter.make_count_file()
    if not return_normalized:
        counter.counts *= counter.std
        counter.counts += counter.mean

    names = get_names_from_counter(counter)
    return (counter.mean, counter.std, counter.counts, names)

def get_names_from_counter(counter):
    if hasattr(counter, 'names') and len(counter.names) > 0:
        names = counter.names
    else:
        names = counter.reader.get_headers()

    return names