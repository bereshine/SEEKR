# -*- coding: utf-8 -*-
"""
session_helper.py

@author: Chris Horning
"""

import os
import pathlib
import uuid
import io

USER_FILE_DIR_ROOT = 'user'

#TODO delete the folder and its contents when the user session times out
def _create_user_directory(user_directory_name):
    """
    create a user folder that will last for the session and hold user created files

    """
    if not os.path.exists(USER_FILE_DIR_ROOT):
        os.mkdir(USER_FILE_DIR_ROOT)

    user_dir_path = _get_user_directory_path(user_directory_name)
    if not os.path.exists(user_dir_path):
        os.mkdir(user_dir_path)


def _get_user_directory_path(user_directory_name):
    return os.path.join(USER_FILE_DIR_ROOT, user_directory_name)

def _get_user_file_path(session, file_identifier, extension=''):
    user_directory_name = session['directory_id']
    user_dir_path = _get_user_directory_path(user_directory_name)
    if len(extension) > 0:
        file_identifier = file_identifier + '.' + extension;
    return os.path.join(user_dir_path, file_identifier)


def create_file(file, session, file_identifier, extension=''):
    file_path = _get_user_file_path(session, file_identifier, extension)

    if isinstance(file, str):
        with open(file_path, 'w') as save_file:
            save_file.write(file)
    else:
        file.save(file_path)



def get_file(session, file_identifier, extension=''):
    file_path = _get_user_file_path(session, file_identifier, extension)
    with open(file_path) as file:
        contents = file.read()

    return contents

def get_file_for_directory_id(directory_id, file_identifier, extension=''):
    mock_session = {'directory_id':directory_id}
    return get_file(mock_session, file_identifier, extension)

def get_directory_id(session):
    return session['directory_id']

def generate_file_identifier():
    return str(uuid.uuid4())


#TODO if a user has a named login, we could change the folder to use that username
def init_user_login(session):
    """
    Initialize any initial data structures for a newly logged in user
    """

    session['directory_id'] = generate_file_identifier()
    _create_user_directory(session['directory_id'])

    session['logged_in'] = True

