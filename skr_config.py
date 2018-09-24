LOGIN_ENABLED = False
CURR_GENCODE_HUMAN_FTP = 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/gencode.v27.lncRNA_transcripts.fa.gz'
CURR_GENCODE_MOUSE_FTP = 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M15/gencode.vM15.lncRNA_transcripts.fa.gz'

MAX_KMER_LENGTH_PRECOMPUTE = 7
LOGGER_LEVEL = 'DEBUG'  #'DEBUG', 'WARN', 'ERROR'

from collections import namedtuple
SavedSet = namedtuple('SavedSet', 'url server_name display_name')
GENCODE_HUMAN = SavedSet(CURR_GENCODE_HUMAN_FTP, 'gencode_human_set', 'All Human (Gencode)')
GENCODE_MOUSE = SavedSet(CURR_GENCODE_MOUSE_FTP, 'gencode_mouse_set', 'All Mouse (Gencode)')

SETTING_USER_SET = 'user_set'
SETTING_COMPARISION_SET = 'comparison_set'

SEQUENCE_NAME_DISPLAY_LENGTH = 20
MAX_VISUAL_SEQ_LENGTH = 200
MAX_VISUAL_KMER_LENGTH = 6