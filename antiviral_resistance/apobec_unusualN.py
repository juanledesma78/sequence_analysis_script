import sys
inFile = sys.argv[1]

import pandas as pd
pd.set_option('display.max_columns', None, 'display.max_rows', None)
data = pd.read_csv(inFile, sep = '\t')
#data = pd.read_csv('./Resistance/SequenceSummaryMOD.tsv', sep = '\t')


#The selection is done using the 20PC sequences as 2PC would have a lot of apobec mutations due to multiple ambiguities 
pc20 = data['Sequence Name'].str.contains('_20PC')
apobec = data['Num Apobec Mutations'] <= 2

#data_pc20 = data[pc20]
#data_filtered_1 = data[pc20].loc[apobec, 'Sequence Name'].str.replace('_20PC','')
data_filtered_1 = data.loc[pc20 & apobec, 'Sequence Name'].str.replace('_20PC','')
# ending _20PC is removed to use just the RS number so a filter will included sequences at 20 and 2PC
# strip removes extra numbers for any reason

# Create a list with the selected sequences that will be used for the final step
selected_seq=[]
for d in range(len(data_filtered_1)):
    
    selected_seq.append(data_filtered_1.values[d])
print(selected_seq)

#selection = data['Sequence Name'].str.contains('701','_20PC')
# It selects only the first condition if exists, no the second
#QC = data[selection]
#print(QC)

# include a new column in the original dataframe just with RS numbers which will be use to select from the list selected_seq
data[['RS', 'Freq']] = data['Sequence Name'].str.split('_', expand = True)
selection = data['RS'].isin(selected_seq)
NO_apobec_sequences = (data[selection]).drop(columns = ['RS', 'Freq'])

NinPRRTIN = (NO_apobec_sequences['PR Other'].str.contains('X')) & (NO_apobec_sequences['RT Other'].str.contains('X')) & (NO_apobec_sequences['IN Other'].str.contains('X'))

#NinPRRTIN = (data[column=(['PR Other','RT Other','IN Other'])].str.contains('X'))
#It should be better to do it independently for each target


final = NO_apobec_sequences.loc[~NinPRRTIN,['Sequence Name','PR Other','RT Other','IN Other']]
#d = data['PR Other'].str.contains('X')
#NinPR = data.loc[d,'PR Other'].str.split(',', expand = True)
print ('These are the sequences with no Ns in the three targets')
print(final['Sequence Name'])

#(data.loc[NinPRRTIN,['Sequence Name','PR Other','RT Other','IN Other', 'Num BDHVN', 'Num Unusual Mutations']]).to_csv('Ns_everywhere.csv', index = False)

#NO_apobec_sequences.to_csv('Sequences_QC_apobec.csv' ,header = True, index = False)
#final.to_csv('Sequences_QC_apobec_Ns.csv' ,header = True, index = False)