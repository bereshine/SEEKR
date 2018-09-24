import urllib.request
import gzip
import os
import shutil

# 1st arg : url link for the file
def download_unzip(saved_set):
    url = saved_set.url
    # create cache directory
    if not os.path.exists("cache"):
        os.mkdir("cache")

    # read file using url
    res = urllib.request.urlopen(url)
    
    # unzip the file into file_content
    file = gzip.open(res, 'rb')
    file_content = file.read()
    file.close()

    # create the unzipped file from file_content
    with open("cache/"+get_unzipped_file_name(saved_set), 'b+w') as unzipped_File:
        unzipped_File.write(file_content)

def get_unzipped_file_name(saved_set):
    url = saved_set.url
    # grab the file name from the url
    index = url.index("gencode.")
    unzipped_file_name = url[index: len(url) - 3]
    return unzipped_file_name

def clear_cache():
    if not os.path.exists('cache'):
        return

    shutil.rmtree('cache')


#download_unzip("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/gencode.v27.lncRNA_transcripts.fa.gz")