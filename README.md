SEEKR Web Portal

Administrator Manual
General Information:
Python 3.5:
Git Repo: sc.unc.edu/chrisrh/skr

This application is dependant upon several libraries for Python 3.5:
pip==9.0.1
Flask==0.12.2
tqdm==4.17.1
numpy==1.13.1

Which are specified in the requirements.txt document

The application’s environment variables are set and initialized by app.py, including --entry-point, which points the server to the main program to execute. 

Other configuration options are specified in skr_config.py, including constants

File Size Requirements:
Files run with and produced by the algorithm can be quite large. The following file size details should give an idea of minimum system requirements that should be met to run SEEKR effectively. 
Gencode human
27908 sequences
519,703 lines
Sequences length:
mean=1027, std=1867, median 636
longest=205012, #10=23112, #300 < 5000
File size – 8MB zipped; 30MB unzipped
Gencode mice
16679 sequences
24MB unzipped
Large test file
234,537 sequences
29,748,112 lines
File size – 560MB zipped; 2GB unzipped
Frequency counts
4^k * float32size * N_sequences
Pearsons
N*N*double_size
8 bytes
Example:
5942MB for gencode pearsons as numpy array
410GB for large test file
~1.8MB for every row of matrix
corr2
N*M*double_size
Cache file sizes of the different kmer lengths
Mice
266,944
1,067,536
4,269,904
17,079,376
68,317,264
273,268,816
1,093,075,024
Human
446,608
1,786,192
7,144,528
28,577,872
114,311,248
457,244,752
1,828,978,768

Running Locally:
Create a local instance of the application on your machine by either cloning the git repository, or downloading it contents

Before spinning the local server up make sure to run precompute_sequence_sets.py. If this script is not run before initializing the application, errors will be encountered when using the standard gencode sets to run the algorithm.

Executing either app.py or seekrServer.py will boot the application locally.


Deployment:
This application is built to run on UNC Cloud Apps’ OpenShift system. 

The first step to deployment is installing Openshift command line tool, instructions for which can be found here: https://help.unc.edu/help/carolina-cloudapps-installing-the-command-line-cli-tools/ 

Once the command line tool is installed, use it to establish an ssh pipeline between the server and the project’s git repository by following the instructions here: https://help.unc.edu/help/carolina-cloudapps-ssh-keys/ 

Once the ssh keys are created and saved, the OpenShift console allows you to build at will and will automatically deploy successful builds to the server. 

Important: precompute_sequence_sets.py must be executed on the server before the application can be run correctly. This can be done from your local machine using the command line tool, ‘os python3 precompute_sequence_sets.py’ or from the web console’s built in terminal. 
