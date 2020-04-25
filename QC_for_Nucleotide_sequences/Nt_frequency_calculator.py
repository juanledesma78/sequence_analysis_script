import sys
inFile = sys.argv[1]
input_file = open(inFile, 'rt')
seq_line=""
headers=[]
sequences=[]

for l in input_file:
    line = l.upper()
    if ">" in line:
        hdr=line.strip('>')
        if len(headers)>0: 
            sequences.append(seq_line)
        seq_line=""
        y=hdr.rstrip()
        headers.append(y)           
    else:
        y=line.rstrip()
        seq_line+=y 
sequences.append(seq_line)

fl_name = inFile.replace('.fas','_Nucelotides_composition.txt')
QC_fl_name = inFile.replace('.fas','_NsFreq_QC_NOT_passed.txt')
summary_file = open(fl_name,'wt')
QC_file = open(QC_fl_name,'wt')
QC_file.write('The frequency of Ns in the following sequences is higher 20% so were excluded from the analysis:'+('\n')*2)

IUPAC_nts = ['A', 'C', 'G', 'T', 'U', 'W', 'S', 'M', 'K', 'R', 'Y' , 'B', 'D', 'H', 'V', 'N']

for i in range(len(headers)):
    summary_title_col= headers[i] +'\t'+'n'+'\t' + 'Frequency'+ '\n'
    summary_file.write(summary_title_col)
    
    for x in range(len(IUPAC_nts)):
        
        if IUPAC_nts[x] in sequences[i]:
            ungapped_sequence = sequences[i].replace('-','') #in case the sequence has gaps (-)
            genome_length = len(ungapped_sequence)
            nt= sequences[i].count(IUPAC_nts[x])
            nt_PC = (nt*100)/genome_length
            nt_info = IUPAC_nts[x] + '\t'+ str(nt) + '\t'+ str("%.2f" % round(nt_PC,2)) +'%' + '\n'
            summary_file.write(nt_info)
            
            if IUPAC_nts[x] =='N':
                if nt_PC > 20: 
                    QC_file.write(headers[i] + '\t'+ str("%.2f" % round(nt_PC,2)) +'%' + '\n')
            
    summary_file.write('lentgh' + '\t'+ str(genome_length) + '\n'+ '\n')
    
summary_file.close()
QC_file.close()

