import os
import mod_wsgi.server

"""
app.py is the file that Apache's mod_wsgi is looking for to setup the environment for the WSGI server.
It specifies the entry point seekrServer.py for the WSGI framework to find the Flask application which
should be named "application". See https://www.python.org/dev/peps/pep-0333/ for the WSGI standard

Notes
-----
Most of the default timeouts were increased to accommodate long requests 


"""

mod_wsgi.server.start(
  '--log-to-terminal',
  '--port', '8080',
  '--trust-proxy-header', 'X-Forwarded-For',
  '--trust-proxy-header', 'X-Forwarded-Port',
  '--trust-proxy-header', 'X-Forwarded-Proto',
  '--processes', os.environ.get('MOD_WSGI_PROCESSES', '1'),
  '--threads', os.environ.get('MOD_WSGI_THREADS', '10'),
  '--application-type', 'module',
  '--entry-point', 'seekrServer',
  '--limit-request-body', '104857600',
  '--request-timeout', '600',
  '--socket-timeout', '600',
  '--header-timeout', '120',
  '--header-max-timeout', '600',
  '--body-timeout', '300',
  '--body-max-timeout', '600',
  '--queue-timeout', '300',
  '--deadlock-timeout', '300',
  '--startup-timeout', '120',
  '--connect-timeout', '300',
  '--queue-timeout', '200',
  '--graceful-timeout', '90',
  '--verbose-debugging'
)
