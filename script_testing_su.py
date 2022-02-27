import argparse
import sequence_utilities as su
import glob
import os

parser = argparse.ArgumentParser(description="Select fasta sequences from a list")
parser.add_argument('-fasta', required=True)
parser.add_argument('-txt', required=True)
args = parser.parse_args()
fasta = args.fasta
txt =args.txt

#su.fasta_selection_by_id(fasta, txt)
#su.fasta_selection_by_strings(fasta)
#su.extract_single_fasta_from_alignment(fasta, False)


#current_directory = os.getcwd()
path_to_INITIO = "testing_data/INITIO/"
for fasta in glob.glob(os.path.join(path_to_INITIO,'*.fas')):
    su.rename_fasta_quasibam_files(fasta)