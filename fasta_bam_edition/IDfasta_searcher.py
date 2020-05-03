import sys
input1 = sys.argv[1]
input2 = sys.argv[2]
FASTA_file = open(input1, 'rt')
ID_file = open(input2,'rt')
seq_line=""
headers=[]
sequences=[]
for line in FASTA_file:
    if ">" in line:
        hdr=((line).strip('>')).rstrip()
        if len(headers)>0: 
            sequences.append(seq_line)
        seq_line=""
        headers.append(hdr)           
    else:
        y=(line.rstrip()).replace('-','')
        seq_line+=y.upper() 
sequences.append(seq_line)
Seq_dictionary = dict(zip(headers, sequences))
fastaID_file = input1.replace('.fas','_IDselection.fas')
outfile = open(fastaID_file,'wt')
for id in ID_file:
    i=id.rstrip()
    if i=='':
        pass
    else:
        for hd,sq in Seq_dictionary.items():
            if i in hd:
                outfile.write('>' + str(hd) +'\n' + str(sq) +'\n')
FASTA_file.close()
ID_file.close()
outfile.close()