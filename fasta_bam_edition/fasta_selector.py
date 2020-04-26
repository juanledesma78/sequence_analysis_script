file = input ('Specify the path/file to analyse:\n')
input_file = open(file, 'rt')
sqs=""
headers=[]
sequences=[]

# Fist line in the loop would be > and the following the sequence of the fasta sequences.
# For the first line with '>' the code jumps the first conditional as header length is 0, and remove \n and append the header.
# Second line would be sequence, strip \n and save in sqs and keeps concatenate info.
# Next '>' line goes to the first statement and there is already a header (len>0),
# firstly apends the concatenated sqs corresponding to the sequences variable, set sqs to zero and appends second header.
# Keeps doing it until the end, but the last concatenation of sequences (sqs) needs to be added extrenally
# otherwise, last lines for the last header would be missed

for line in input_file:

    if ">" in line: 
        if len(headers)>0: 
            sequences.append(sqs)
        sqs=""
        y=line.rstrip()
        headers.append(y)           
    else:
        y=line.rstrip()
        sqs+=y # to get all the lines from the same sequences in a single string
          
sequences.append(sqs) 


# This code creates to two individual fasta files with all the 20PC and 2PC sequences, respectively
selected_20PC=[]
for n in range(len(headers)):
    if '_20PC' in headers[n]:
        y = headers[n]+'\n'+ sequences[n]+'\n'
        selected_20PC.append(y)
        
x="".join(selected_20PC)
name=file.replace(".fas","_20PC.fas") #using strip the first character is missing for any reason , so better use replace
output_file= open(name,'w')
output_file.write(x)
output_file.close()

selected_2PC=[]
for n in range(len(headers)):
    if '_2PC' in headers[n]:
        y = headers[n]+'\n'+ sequences[n]+'\n'
        selected_2PC.append(y)
        
x="".join(selected_2PC)
name=file.replace(".fas","_2PC.fas")
output_file= open(name,'w')
output_file.write(x)
output_file.close()


# This code creates a single fasta file for each sequence in the alignment. 
for n in range(len(headers)):
   
        info = headers[n]+'\n'+ sequences[n]+'\n'
        file_name = (headers[n]).strip('>')+'.fas'
        output_file= open(file_name,'w')
        output_file.write(info)
        output_file.close()

input_file.close()

