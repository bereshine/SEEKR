"""
Web Server for SEEKR Web Portal using Flask

@author: Chris Horning, Shuo Wang, Qidi Chen

"""
import email.utils
import itertools
import os
import time
import zipfile
from io import BytesIO
from logging.handlers import RotatingFileHandler

import numpy as np
from flask import Flask
from flask import jsonify
from flask import redirect
from flask import render_template
from flask import request
from flask import session
from scipy.stats import stats

import cluster_vis
import session_helper
import skr_config
from SeekrServerError import SeekrServerError
from pearson import pearson
from precompute_sequence_sets import initialize_cache
from seekrLauncher import _run_seekr_algorithm
from seekrLauncher import fixup_counts
from seekrLauncher import get_kmers_csv
from seekrLauncher import get_pearsons_csv

"""
seekrServer.py contains the Web Services the application provides using the Flask framework.

"""

# create app instance
application = Flask(__name__)
application.config['MAX_CONTENT_LENGTH'] = 50 * 1024 * 1024
application.config['SECRET_KEY'] = os.urandom(12)
if not application.debug:
    file_handler = RotatingFileHandler('seekr_server.log')
    file_handler.setLevel(skr_config.LOGGER_LEVEL)
    application.logger.addHandler(file_handler)
application.logger.setLevel(skr_config.LOGGER_LEVEL)

# route handling
@application.route('/')
def admin():
    return redirect('/login')

# login page
@application.route('/login',methods=['POST','GET'])
def login():
    if skr_config.LOGIN_ENABLED is False:
        session_helper.init_user_login(session)
        return redirect('/home')

    if request.method == 'GET':
        if session.get('logged_in')==True:
            return redirect('/home')
        return render_template('login.html')
    if request.method == 'POST':
        if request.form['password'] == '123' and request.form['username'] == '123' \
                or request.form['password'] == 'qidi' and request.form['username'] == 'skr' \
                or request.form['password'] == 'david' and request.form['username'] == 'skr' \
                or request.form['password'] == 'chris' and request.form['username'] == 'skr':
            session_helper.init_user_login(session)
            return redirect('/home')
        else:
            return redirect('/login')

# home page
@application.route('/home',methods=['POST','GET'])
def home():
    if request.method == 'GET':
        if session.get('logged_in') == True:
            return render_template('home.html')
        else:
            return redirect('/login')
    if request.method == 'POST':
        print('checkbox selected')
        print (request.form.getlist('check'))
        print('dropdown selected')
        print (request.form.getlist('drop'))

        return render_template('home.html')

def return_file(file_contents, pearsons):
    t1 = time.perf_counter()
    #Create an in-memory zip file with GZIP to return
    bytes_io = BytesIO()
    zip = zipfile.ZipFile(bytes_io, mode='w', compression=zipfile.ZIP_DEFLATED)
    zip.writestr('counts.csv', str.encode(file_contents, 'utf-8'))
    zip.writestr('pearsons.npy', pearsons)
    zip.close()
    bytes_io.seek(0)
    zipped_bytes = bytes_io.read()
    last_modified = email.utils.formatdate(time.time(), usegmt=True)
    headers = {'Content-Type': 'application/zip',
               'Content-Disposition': 'attachment;filename = seekr.zip',
               'Content-Length' : str(len(zipped_bytes)),
               'Last-Modified' : last_modified
               }
    t2 = time.perf_counter()
    application.logger.debug('Zipping file of length ' + str(len(zipped_bytes))  + ' took %.3f seconds' % (t2 - t1))
    return (zipped_bytes, headers)

# routing function for processing upload actions
@application.route('/jobs', methods=['POST'])
def process_jobs():
    try:
        if skr_config.LOGIN_ENABLED and session.get('logged_in') != True:
            return redirect('/login')

        parameters = build_seekr_parameters(request)
        t1 = time.perf_counter()
        counts, names, comparison_counts, comparison_names, counter = _run_seekr_algorithm(parameters=parameters)
        t2 = time.perf_counter()
        application.logger.debug('Running the algorithm took %.3f seconds' % (t2 - t1))

        if len(names) <= skr_config.MAX_VISUAL_SEQ_LENGTH and len(comparison_names) <= skr_config.MAX_VISUAL_SEQ_LENGTH and parameters['kmer_length'] <= skr_config.MAX_VISUAL_KMER_LENGTH:

            fixup_counts_warnings = fixup_counts(counts, counter)
            if comparison_counts is None:
                comparison_counts = counts
                comparison_names = names
            else:
                fixup_comparision_warnings = fixup_counts(comparison_counts, counter)

            #reorder according to hierarchical cluster example
            Z = cluster_vis.cluster_kmers(counts)
            ordering = cluster_vis.get_ordering(Z)
            ordered_counts = counts[ordering,:]
            ordering_int_list = ordering.astype(int).tolist()
            ordered_names = [names[i] for i in ordering_int_list]

            comparison_Z = cluster_vis.cluster_kmers(comparison_counts)
            comparison_ordering = cluster_vis.get_ordering(comparison_Z)
            comparison_ordered_counts = comparison_counts[comparison_ordering, :]
            comparison_ordering_int_list = comparison_ordering.astype(int).tolist()
            comparison_ordered_names = [comparison_names[i] for i in comparison_ordering_int_list]

            pearsons = pearson(counts, comparison_counts)

            # shorten length of names returned down to 20 characters
            new_names = []
            for s in names:
                if len(s) > skr_config.SEQUENCE_NAME_DISPLAY_LENGTH:
                    new_names.append(s[:skr_config.SEQUENCE_NAME_DISPLAY_LENGTH])
                else:
                    new_names.append(s)

            names = new_names

            new_names = []
            for s in comparison_names:
                if len(s) > skr_config.SEQUENCE_NAME_DISPLAY_LENGTH:
                    new_names.append(s[:skr_config.SEQUENCE_NAME_DISPLAY_LENGTH])
                else:
                    new_names.append(s)

            comparison_names = new_names

            x = ['A', 'G', 'T', 'C']
            kmer = [p for p in itertools.product(x, repeat=parameters['kmer_length'])]

            count = 0
            for i in kmer:
                kmer[count] = ''.join(i)
                count = count + 1

            norm_npm = counts
            flat_npm = norm_npm.flatten()
            scale_npm = norm_npm.flatten()
            mean = np.mean(scale_npm)
            z_npm = stats.zscore(flat_npm)
            count = 0
            for i in z_npm:
                if i >= 2:
                    scale_npm[count] = mean
                elif i < -1:
                    scale_npm[count] = mean
                count = count + 1
            clean_counts = np.reshape(scale_npm, np.shape(norm_npm))

            pearsons = str(pearsons.tolist())
            counts = str(counts.tolist())
            clean_counts = str(clean_counts.tolist())

            return jsonify({'user_names': names, 'comparison_names': comparison_names,
                            'kmer_bins': kmer, 'pearson_matrix': pearsons, 'kmer_matrix': counts,
                            'kmer_matrix_clean': clean_counts, 'user_cluster': ordering_int_list,
                            'comparison_cluster': comparison_ordering_int_list, 'user_warnings': fixup_counts_warnings,
                            'comparison_warnings': fixup_comparision_warnings
                            })

        else:
            return jsonify({'visual_flag': True})



    except Exception as e:
        application.logger.exception('Error in /jobs')
        return jsonify({'error': "Server_Error: 500"})

# routing function for generating kmers file
@application.route('/files/kmers', methods=['POST'])
def process_kmer_job():
    try:
        if skr_config.LOGIN_ENABLED and session.get('logged_in') != True:
            return redirect('/login')

        parameters = build_seekr_parameters(request)

        application.logger.debug(parameters)

        t1 = time.perf_counter()
        counts, names, comparison_counts, comparison_names, counter = _run_seekr_algorithm(parameters=parameters)
        t2 = time.perf_counter()
        application.logger.debug('Running the algorithm took %.3f seconds' % (t2 - t1))

        fixup_counts_warnings = fixup_counts(counts, counter)
        if comparison_counts is None:
            comparison_counts = counts
            comparison_names = names
        else:
            fixup_comparision_warnings = fixup_counts(comparison_counts, counter)


        x = ['A', 'G', 'T', 'C']
        kmer = [p for p in itertools.product(x, repeat=parameters['kmer_length'])]

        csv_string = get_kmers_csv(counts=counts, names=names, kmers=kmer)

        last_modified = email.utils.formatdate(time.time(), usegmt=True)
        headers = {'Content-Type': 'application/csv',
                   'Content-Disposition': 'attachment;filename = seekr.csv',
                   'Content-Length': str(len(csv_string)),
                   'Last-Modified': last_modified
                   }
        return (csv_string, headers)

    except Exception as e:
        application.logger.exception('Error in /files/kmers')
        #TODO change error from json
        return jsonify({'error': "Server_Error: 500"})

# routing function for generating pearsons file
@application.route('/files/pearsons', methods=['POST'])
def process_pearsons_job():
    try:
        if skr_config.LOGIN_ENABLED and session.get('logged_in') != True:
            return redirect('/login')

        parameters = build_seekr_parameters(request)

        application.logger.debug('CURRENT METHOD: process_pearsons_job')

        t1 = time.perf_counter()
        counts, names, comparison_counts, comparison_names, counter = _run_seekr_algorithm(parameters=parameters)
        t2 = time.perf_counter()
        application.logger.debug('Running the algorithm took %.3f seconds' % (t2 - t1))

        fixup_counts_warnings = fixup_counts(counts, counter)
        if comparison_counts is None:
            comparison_counts = counts
            comparison_names = names
        else:
            fixup_comparision_warnings = fixup_counts(comparison_counts, counter)

        pearsons = pearson(counts, comparison_counts)
        csv_string = get_pearsons_csv(names, pearsons, comparison_names)

        last_modified = email.utils.formatdate(time.time(), usegmt=True)
        headers = {'Content-Type': 'application/csv',
                   'Content-Disposition': 'attachment;filename = pearsons.csv',
                   'Content-Length': str(len(csv_string)),
                   'Last-Modified': last_modified
                   }
        return (csv_string, headers)

    except Exception as e:
        application.logger.exception('Error in /files/pearsons')
        #TODO change error from json
        return jsonify({'error': "Server_Error: 500"})

def build_seekr_parameters(request):
    # TODO provide reasonable defaults

    parameters = dict()

    if(request.get_json()):
        json_parameters = request.get_json()

        if ('user_set_id' in json_parameters):
            parameters['user_set_files'] = json_parameters['user_set_id']
        if ('comparison_set_id' in json_parameters and json_parameters['comparison_set_id'] is not None and
                    json_parameters[
                        'comparison_set_id'] != ''):
            parameters['comparison_set_files'] = json_parameters['comparison_set_id']
        if 'comparison_set' in json_parameters:
            parameters['comparison_set'] = str(json_parameters['comparison_set'])
        if 'kmer_length' in json_parameters:
            parameters['kmer_length'] = int(json_parameters['kmer_length'])
        if 'normal_set' in json_parameters:
            parameters['normal_set'] = str(json_parameters['normal_set'])
        parameters['directory_id'] = session_helper.get_directory_id(session)
        if parameters['directory_id'] is None or len(parameters['directory_id']) <= 0:
            raise SeekrServerError('User directory not found for this session')

    else:

        if ('user_set_id' in request.form):
            parameters['user_set_files'] = request.form['user_set_id']
        if ('comparison_set_id' in request.form and request.form['comparison_set_id'] is not None and
                    request.form[
                        'comparison_set_id'] != ''):
            parameters['comparison_set_files'] = request.form['comparison_set_id']
        if 'comparison_set' in request.form:
            parameters['comparison_set'] = str(request.form['comparison_set'])
        if 'kmer_length' in request.form:
            parameters['kmer_length'] = int(request.form['kmer_length'])
        if 'normal_set' in request.form:
            parameters['normal_set'] = str(request.form['normal_set'])
        parameters['directory_id'] = session_helper.get_directory_id(session)
        if parameters['directory_id'] is None or len(parameters['directory_id']) <= 0:
            raise SeekrServerError('User directory not found for this session')

    return parameters


def check_sequence_length(file_str, upper_bound):
    count = file_str.count('>')
    print(count)
    file_more_than_200_sequences = (count / 2) > upper_bound
    return file_more_than_200_sequences


@application.route('/files/fasta', methods=['POST'])
def create_fasta():
    """
    Post to upload a new fasta file
    This file will be given a unique identifier

    """
    assert request.method == 'POST'

    if 'file' not in request.files:
        application.logger.debug('Error, no file')
        # TODO error case

    file = request.files['file']

    file_identifier = session_helper.generate_file_identifier()
    session_helper.create_file(file, session, file_identifier, extension='fasta')

    # get uploaded file from session
    get_file = session_helper.get_file(session, file_identifier, 'fasta')

    file_more_than_200_sequences = check_sequence_length(get_file, skr_config.MAX_VISUAL_SEQ_LENGTH)

    json_dict = {
        'file_id': file_identifier,
        'file_more_than_200_sequences': file_more_than_200_sequences
    }

    return jsonify(json_dict)


if __name__ == '__main__':
    application.run()

