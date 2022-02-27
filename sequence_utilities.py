from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import re
import glob

def fasta_selection_by_id(fasta, txt):
    """It takes a given fasta file with several sequences and returns a new 
    fasta file containing the sequences matching the ids given in txt file"""

    fasta_input = os.path.join(os.getcwd(),fasta)
    seq_list = os.path.join(os.getcwd(),txt)

    with open(seq_list,'r') as seqid:
        seqidlist = [line.rstrip() for line in seqid ]
    
    selection = {} # remove duplicates
    for record in SeqIO.parse(fasta_input, 'fasta'): #Bio.SeqIO.FastaIO.FastaIterator object
        if record.id in seqidlist:
            selection[record.id] = SeqRecord(record.seq, record.id, description='')
    output_file = fasta_input.replace('.fas','_selection_IDs.fasta')
    SeqIO.write(selection.values(), output_file,'fasta')


def fasta_selection_by_strings(fasta):
    """It takes a given fasta file and and returns new fasta file(s) containing 
    the sequences which headers match the strings entered in a prompt"""

    search = input("Enter the string(s) to find in the headers separated ONLY by comma: \n")
    searching_list = search.split(',')
    fasta_input = os.path.join(os.getcwd(),fasta)
    for n in range(len(searching_list)):
        selection = {} # remove duplicates
        for record in SeqIO.parse(fasta_input, 'fasta'): 
            if searching_list[n] in record.id:
                selection[record.id] = SeqRecord(record.seq, record.id, description='')
        if selection.values != '':
            output_file = fasta_input.replace('.fas',f'_selection_{searching_list[n]}.fasta')
            SeqIO.write(selection.values(), output_file,'fasta')


def extract_single_fasta_from_alignment(fasta, ungap=True):
    """It takes a given FASTA file and returns as many new FASTA files as 
    sequences contained in the given file. By default, the gaps are removed
    from the sequences in the new files"""

    fasta_input = os.path.join(os.getcwd(),fasta)
    output_path = os.path.dirname(fasta_input)
    for record in SeqIO.parse(fasta_input, 'fasta'): 
        if ungap == True:
            sequence = record.seq.ungap('-')
            fasta_output = SeqRecord(sequence, record.id, description='')
            SeqIO.write(fasta_output, 
                        os.path.join(output_path,f'{record.id}.fasta'),
                        'fasta')
        else:
            fasta_output = SeqRecord(record.seq, record.id, description='')
            SeqIO.write(fasta_output, 
                        os.path.join(output_path,f'{record.id}_gapped.fasta'),
                        'fasta')

def rename_fasta_quasibam_files(fasta_input):
    """INITIO project. It takes the FASTA files generated from a NGS pipeline
    and rename the FASTA files, use the same notation for the FASTA headers and
    rename QUASIBAM files"""

    fasta = SeqIO.read(fasta_input, 'fasta')
    id = (re.sub(r'^\d+_','', fasta.id)).replace('%_','PC')
    fasta.id = id
    fasta.description = ''
    fasta_output = re.sub(r'\d+_[H|RS].+\.fas',f'{id}.fasta', fasta_input)
    SeqIO.write(fasta, fasta_output, 'fasta')
    
    tabular_input = fasta_input.replace('.20.fas', '.tabular')
    #bam_input = tabular_input.replace('.tabular','.bam')
    if tabular_input: #and bam_input:
        if '.20.fas' in fasta_input and '_20PC.fasta' in fasta_output:
            tabular_output = fasta_output.replace('_20PC.fasta', '.tabular')
            os.rename(tabular_input, tabular_output)
            #bam_output = tabular_output.replace('.tabular','.bam')
            #os.rename(bam_input, bam_output)
    os.remove(fasta_input)


def NcontentQC():
    pass
# #NcontentQC_FINALv1.py
# import sys
# import re
# inFile = sys.argv[1]
# #cutoff = sys.argv[2] #in case the cutoff is taken as argument at the time of running the script
# input_file = open(inFile, 'rt')
# 
# """reading fasta file"""
# seq_line=""
# h=[]
# s=[]
# for line in input_file:
#     if ">" in line:
#         hdr=(line.strip('>')).rstrip() 
#         if len(h)>0: 
#             s.append(seq_line)
#         seq_line=""
#         h.append(hdr)           
#     else:
#         y=(line.rstrip()).replace('-','')
#         seq_line+=y.upper() 
# s.append(seq_line)
# 
# """Selecting the correct sequences and filtering out negative controls, no specified IDs"""
# idx_allowed_list = [i for i, item in enumerate(h) if re.search('(RS\d{8,10}|H\d{9,11}|NC_045512.2)', item)]
# headers = [h[n] for n in idx_allowed_list]
# sequences = [s[n] for n in idx_allowed_list]
# 
# 
# """QC and good quality sequences saved in a new fasta file"""
# QC_name = inFile.replace('.fas','_Nucleotides_composition_QC_v1.txt')
# QC_summary = open(QC_name,'wt')
# FASTA_name = inFile.replace('.fas','_QC_passed_v1.fas')# regular expression to take fas.fasta.txt ?
# FASTA_selected = open(FASTA_name,'wt')
# comments = '\n'+'N% report:'+'\n'+'The frequency of undetermined nucleotides (Ns) in the following sequences is higher than 20% so they will be excluded from further analysis'+('\n')*2
# IUPAC_nts = ['A', 'C', 'G', 'T', 'U', 'W', 'S', 'M', 'K', 'R', 'Y' , 'B', 'D', 'H', 'V', 'N']
# 
# seq_QC_passed = ""
# for i in range(len(headers)):
#     title_columns= headers[i] +'\t'+'n'+'\t' + 'Frequency'+ '\n'
#     QC_summary.write(title_columns)
#     try:
#         for x in range(len(IUPAC_nts)):
#             if IUPAC_nts[x] in sequences[i]:
#                 genome_length = len(sequences[i])
#                 nt= sequences[i].count(IUPAC_nts[x])
#                 nt_PC = (nt*100.00)/genome_length
#                 nt_info = IUPAC_nts[x] + '\t'+ str(nt) + '\t'+ str("%.2f" % round(nt_PC,2)) +'%' + '\n'
#                 QC_summary.write(nt_info)
#                 
#                 if IUPAC_nts[x] =='N':
#                     #if nt_PC <= int(cutoff):
#                     if nt_PC <= 20:
#                         h = '>' + headers[i] + '\n'+ sequences[i] + '\n'
#                         seq_QC_passed += h
#                     if nt_PC > 20:
#                         qc=headers[i]  + '\t'+ str("%.2f" % round(nt_PC,2)) +'% Ns in sequence' + '\n'
#                         comments +=qc
#         if 'N' not in sequences[i]:
#             h = '>'+ headers[i] + '\n'+ sequences[i] + '\n'
#             seq_QC_passed += h 
#         QC_summary.write('lentgh' + '\t'+ str(genome_length) + '\n'+ '\n')
#     except:
#          QC_summary.write('WARNING: there is no nucleotide data for this sequence, the sequence field from teh fasta file was empty'+'\n'*2)
# FASTA_selected.write(seq_QC_passed)
# QC_summary.write(comments)
# QC_summary.close()
# FASTA_selected.close()
# 
# input_file.close()

    
