import sys
inFile = sys.argv[1]
#cutoff = sys.argv[2] #in case the cutoff is taken as argument at the time of running the script
input_file = open(inFile, 'rt')
QC_name = inFile.replace('.fas','_Nucleotides_composition_QC.txt')
FASTA_name = inFile.replace('.fas','_QC_passed.fas')
QC_summary = open(QC_name,'wt')
FASTA_selected = open(FASTA_name,'wt')

seq_line=""
headers=[]
sequences=[]
comments = '\n'+'WARNING!'+'\n'+'WARNING!'+'\n'+'The frequency of undetermined nucleotides (Ns) in the following sequences is higher than 20% so they will be excluded from further analysis'+('\n')*2
for line in input_file:
    if ">" in line:
        hdr=(line.strip('>')).rstrip()
        if len(headers)>0: 
            sequences.append(seq_line)
        seq_line=""
        headers.append(hdr)           
    else:
        y=(line.rstrip()).replace('-','')
        seq_line+=y.upper() 
sequences.append(seq_line)


IUPAC_nts = ['A', 'C', 'G', 'T', 'U', 'W', 'S', 'M', 'K', 'R', 'Y' , 'B', 'D', 'H', 'V', 'N']
seq_QC_passed = ""
for i in range(len(headers)):
    title_columns= headers[i] +'\t'+'n'+'\t' + 'Frequency'+ '\n'
    QC_summary.write(title_columns)
    for x in range(len(IUPAC_nts)):
        if IUPAC_nts[x] in sequences[i]:
            genome_length = len(sequences[i])
            nt= sequences[i].count(IUPAC_nts[x])
            nt_PC = (nt*100)/genome_length
            nt_info = IUPAC_nts[x] + '\t'+ str(nt) + '\t'+ str("%.2f" % round(nt_PC,2)) +'%' + '\n'
            QC_summary.write(nt_info)
            if IUPAC_nts[x] =='N':
                #if nt_PC > int(cutoff):
                if nt_PC <= 20:
                    h = '>'+ headers[i] + '\n'+ sequences[i] + '\n'
                    seq_QC_passed += h
                if nt_PC > 20:
                    qc=headers[i] + '\t'+ str("%.2f" % round(nt_PC,2)) +'%' + '\n'
                    comments +=qc
                
    QC_summary.write('lentgh' + '\t'+ str(genome_length) + '\n'+ '\n')
FASTA_selected.write(seq_QC_passed)
QC_summary.write(comments)
QC_summary.close()
FASTA_selected.close()
input_file.close()
