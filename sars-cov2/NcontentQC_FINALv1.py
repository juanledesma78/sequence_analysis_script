import sys
import re
inFile = sys.argv[1]
#cutoff = sys.argv[2] #in case the cutoff is taken as argument at the time of running the script
input_file = open(inFile, 'rt')

"""reading fasta file"""
seq_line=""
h=[]
s=[]
for line in input_file:
    if ">" in line:
        hdr=(line.strip('>')).rstrip() 
        if len(h)>0: 
            s.append(seq_line)
        seq_line=""
        h.append(hdr)           
    else:
        y=(line.rstrip()).replace('-','')
        seq_line+=y.upper() 
s.append(seq_line)

"""Selecting the correct sequences and filtering out negative controls, no specified IDs"""
idx_allowed_list = [i for i, item in enumerate(h) if re.search('(RS\d{8,10}|H\d{9,11}|NC_045512.2)', item)]
headers = [h[n] for n in idx_allowed_list]
sequences = [s[n] for n in idx_allowed_list]


"""QC and good quality sequences saved in a new fasta file"""
QC_name = inFile.replace('.fas','_Nucleotides_composition_QC_v1.txt')
QC_summary = open(QC_name,'wt')
FASTA_name = inFile.replace('.fas','_QC_passed_v1.fas')# regular expression to take fas.fasta.txt ?
FASTA_selected = open(FASTA_name,'wt')
comments = '\n'+'N% report:'+'\n'+'The frequency of undetermined nucleotides (Ns) in the following sequences is higher than 20% so they will be excluded from further analysis'+('\n')*2
IUPAC_nts = ['A', 'C', 'G', 'T', 'U', 'W', 'S', 'M', 'K', 'R', 'Y' , 'B', 'D', 'H', 'V', 'N']

seq_QC_passed = ""
for i in range(len(headers)):
    title_columns= headers[i] +'\t'+'n'+'\t' + 'Frequency'+ '\n'
    QC_summary.write(title_columns)
    try:
        for x in range(len(IUPAC_nts)):
            if IUPAC_nts[x] in sequences[i]:
                genome_length = len(sequences[i])
                nt= sequences[i].count(IUPAC_nts[x])
                nt_PC = (nt*100.00)/genome_length
                nt_info = IUPAC_nts[x] + '\t'+ str(nt) + '\t'+ str("%.2f" % round(nt_PC,2)) +'%' + '\n'
                QC_summary.write(nt_info)
                
                if IUPAC_nts[x] =='N':
                    #if nt_PC <= int(cutoff):
                    if nt_PC <= 20:
                        h = '>' + headers[i] + '\n'+ sequences[i] + '\n'
                        seq_QC_passed += h
                    if nt_PC > 20:
                        qc=headers[i]  + '\t'+ str("%.2f" % round(nt_PC,2)) +'% Ns in sequence' + '\n'
                        comments +=qc
        if 'N' not in sequences[i]:
            h = '>'+ headers[i] + '\n'+ sequences[i] + '\n'
            seq_QC_passed += h 
        QC_summary.write('lentgh' + '\t'+ str(genome_length) + '\n'+ '\n')
    except:
         QC_summary.write('WARNING: there is no nucleotide data for this sequence, the sequence field from teh fasta file was empty'+'\n'*2)
FASTA_selected.write(seq_QC_passed)
QC_summary.write(comments)
QC_summary.close()
FASTA_selected.close()

input_file.close()
