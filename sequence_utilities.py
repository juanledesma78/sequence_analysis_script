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


def NcontentQC(fasta):
    fasta_input = os.path.join(os.getcwd(),fasta)
    sequence_input = SeqIO.parse(fasta_input, 'fasta')
    QCPassed = []
    QCNotPassed = []
    reportQC = 'Sequence Id\t N Count\t %\n'
    for record in sequence_input:
        genome_length = len(record.seq)
        Ncount = record.seq.count('N')
        NcountPercentage = (Ncount*100.00)/genome_length
        if NcountPercentage <= 20:
            QCPassed.append(record)
            reportQC += f'{record.id}\t{Ncount}\t{round(NcountPercentage, 2)}\n'
        else:
            QCNotPassed.append(record)
            reportQC += f'{record.id}\t{Ncount}\t{round(NcountPercentage, 2)}\n'

    report_name = re.sub(r'\.(fas|fasta)','_Ncontent_report.txt',fasta_input)
    with open(report_name,'w') as log:
        log.write(reportQC)
    seqsQCpassed = re.sub(r'\.(fas|fasta)','_QC_N%_passed.fasta',fasta_input)
    seqsQCNotpassed = re.sub(r'\.(fas|fasta)','_QC_N%_NOT_passed.fasta',fasta_input)
    SeqIO.write( QCPassed, seqsQCpassed, 'fasta')
    SeqIO.write( QCNotPassed, seqsQCNotpassed, 'fasta')
