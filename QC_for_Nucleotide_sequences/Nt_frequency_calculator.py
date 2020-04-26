import sys
inFile = sys.argv[1]
#cutoff = sys.argv[2] #in case the cutoff is taken as argument at the time of running the script
input_file = open(inFile, 'rt')
output_name = inFile.replace('.fas','_Nucleotides_composition_QC.txt')
QC_summary= open(output_name,'wt')

seq_line=""
headers=[]
sequences=[]
comments = '\n'+'WARNING!'+'\n'+'The frequency of undetermined nucleotides (Ns) in the following sequences is higher than 20% so they will be excluded from further analysis:'+('\n')*2
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
seqs_QC_not_passed = [] #kept to update the script with a tool to select the fasta sequences passing the QC
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
                if nt_PC > 20:
                    qc=headers[i] + '\t'+ str("%.2f" % round(nt_PC,2)) +'%' + '\n'
                    comments +=qc
                    seqs_QC_not_passed.append(headers[i])
                
    QC_summary.write('lentgh' + '\t'+ str(genome_length) + '\n'+ '\n')
QC_summary.write(comments)
QC_summary.close()
input_file.close()
