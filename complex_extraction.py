from __future__ import unicode_literals, print_function
import sys,codecs
import os,csv,re,json


json_file='./Glygen/OGER/uniprot_human_sprot_two_dic.json'
with open(json_file) as jfile:
    uniprot_kb_dic=json.load(jfile)

uniprot_kb_protein_dic=uniprot_kb_dic[0]
uniprot_kb_gene_dic=uniprot_kb_dic[1]

complex_dic={}
#for subunit match
subunit_pattern_1='(\d+\s?kDa)?\s?([Rr]egulatory|[Ee]nzymatic|[Cc]atalytic).*subunit'
subunit_pattern_2='(\d+\s?kDa)?\s?subunit'

#for chain match
chain_pattern_1='(\d+\s?kDa)?\s?([Rr]egulatory|[Ee]nzymatic|[Cc]atalytic).*\schain'
chain_pattern_2='(\d+\s?kDa)?\schain'
#for kDa match
kDa_pattern='\d+\s?kDa'

#for the information after subunit
greek_list=['[Aa]lpha','[Bb]eta','[Gg]amma','[Dd]elta','[Ss]igma']
greek_list_short=[gi[1] for gi in greek_list]
after_subunit_greek='('+'|'.join(greek_list)+')'
after_subunit_greek_digit='('+'|'.join(greek_list)+')-?\s?(\d+)([A-Z]?)'

after_subunit_digit='\d+'
after_subunit_cap_letter='[A-Z]{1,3}'

#comple the patterns
complie_subunit_pattern_1=re.compile(subunit_pattern_1)
complie_subunit_pattern_2=re.compile(subunit_pattern_2)

complie_chain_pattern_1=re.compile(chain_pattern_1)
complie_chain_pattern_2=re.compile(chain_pattern_2)

complie_kDa_pattern=re.compile(kDa_pattern)

complie_after_subunit_greek=re.compile(after_subunit_greek)
complie_after_subunit_greek_digit=re.compile(after_subunit_greek_digit)

complie_after_subunit_single_digit=re.compile(after_subunit_digit)
complie_after_subunit_single_letter=re.compile(after_subunit_cap_letter)

protein_kda_dic={}
protein_subunit_info_dic={}
for ki in uniprot_kb_protein_dic.keys():
    for protein_name in uniprot_kb_protein_dic[ki]:

        sr_subunit_1=complie_subunit_pattern_1.search(protein_name)
        sr_subunit_2=complie_subunit_pattern_2.search(protein_name)

        sr_chain_1=complie_chain_pattern_1.search(protein_name)
        sr_chain_2=complie_chain_pattern_2.search(protein_name)
        sr_kda=complie_kDa_pattern.search(protein_name)

        #find the kda digit to remove relevant the gene name
        #[True/False, kDa digit]
        protein_kda_dic[ki]=[False,'']
        #[greek,greek+digit , digit, single letter]
        if ki not in protein_subunit_info_dic:
            protein_subunit_info_dic[ki]=['','','','']
        if sr_kda:

            kda_digit=protein_name[sr_kda.span()[0]:sr_kda.span()[1]-3]
            kda_digit=kda_digit.strip()
            protein_kda_dic[ki]=[True,kda_digit]

        if sr_subunit_1 or sr_chain_1:
            if sr_subunit_1 and sr_chain_1:
                if sr_subunit_1.span()[1]>sr_chain_1.span()[1]:
                    sr_1=sr_chain_1
                else:
                    sr_1=sr_subunit_1
            elif sr_subunit_1:
                sr_1=sr_subunit_1
            elif sr_chain_1:
                sr_1=sr_chain_1

            if ki in complex_dic:
                complex_dic[ki].append([ki,protein_name,protein_name[:sr_1.span()[0]]])
            else:
                complex_dic[ki]=[[ki,protein_name,protein_name[:sr_1.span()[0]]]]
            #print('ID:',ki)
            #print('Original protein name: ',protein_name)
            #print(sr_subunit_1.span())
            #print('Complex name: ',protein_name[:sr_subunit_1.span()[0]])
            after_subunit_str=protein_name[sr_1.span()[1]+1:]
            after_subunit_str=after_subunit_str.strip()
            sr_greek=complie_after_subunit_greek.search(after_subunit_str)
            sr_greek_digit=complie_after_subunit_greek_digit.search(after_subunit_str)
            
            if sr_greek:
                protein_subunit_info_dic[ki][0]=sr_greek.group(0)
            if sr_greek_digit:
                protein_subunit_info_dic[ki][1]=sr_greek_digit.group(1)[0].upper()+sr_greek_digit.group(2)
            if after_subunit_str.isdigit():
                protein_subunit_info_dic[ki][2]=after_subunit_str
            if len(after_subunit_str)==1 and not after_subunit_str.isdigit():
                protein_subunit_info_dic[ki][3]=after_subunit_str

        elif sr_subunit_2 or sr_chain_2:
            if sr_subunit_2 and sr_chain_2:
                if sr_subunit_2.span()[1]>sr_chain_2.span()[1]:
                    sr_2=sr_chain_2
                else:
                    sr_2=sr_subunit_2
            elif sr_subunit_2:
                sr_2=sr_subunit_2
            elif sr_chain_2:
                sr_2=sr_chain_2

            if ki in complex_dic:
                complex_dic[ki].append([ki,protein_name,protein_name[:sr_2.span()[0]]])
            else:
                complex_dic[ki]=[[ki,protein_name,protein_name[:sr_2.span()[0]]]]
            #print('ID:',ki)
            #print('Original protein name: ',protein_name)
            #print(sr_subunit_2.span())
            #print('Complex name: ',protein_name[:sr_subunit_2.span()[0]])
            after_subunit_str=protein_name[sr_2.span()[1]+1:]
            after_subunit_str=after_subunit_str.strip()
            sr_greek=complie_after_subunit_greek.search(after_subunit_str)
            sr_greek_digit=complie_after_subunit_greek_digit.search(after_subunit_str)
            '''
            if ki=='Q16558':
                print(after_subunit_str)
                if sr_greek_digit:
                    print(sr_greek_digit)
                    print('Group 1:',sr_greek_digit.group(1))
                    print('Group 2:',sr_greek_digit.group(2))
                    print(sr_greek_digit.group(1)[0].upper()+sr_greek_digit.group(2))
            '''
            if sr_greek:
                protein_subunit_info_dic[ki][0]=sr_greek.group(0)
            if sr_greek_digit:
                protein_subunit_info_dic[ki][1]=sr_greek_digit.group(1)[0].upper()+sr_greek_digit.group(2)
            if after_subunit_str.isdigit():
                protein_subunit_info_dic[ki][2]=after_subunit_str
            if len(after_subunit_str)==1 and not after_subunit_str.isdigit():
                protein_subunit_info_dic[ki][3]=after_subunit_str


#for test purpose
uniprot_id='P0DN86'
if uniprot_id in uniprot_kb_protein_dic:
    print(uniprot_kb_protein_dic[uniprot_id])
if uniprot_id in protein_subunit_info_dic:
    print(protein_subunit_info_dic[uniprot_id])
if uniprot_id in uniprot_kb_gene_dic:
    print(uniprot_kb_gene_dic[uniprot_id])
#for the information after subunit
gene_greek_pattern='('+'|'.join(greek_list_short)+')\d*$'

gene_single_digit_pattern='\d+$'
gene_single_letter_pattern='[A-Z]$'
for ki in uniprot_kb_gene_dic.keys():
    for gene_name in uniprot_kb_gene_dic[ki]:
        if ki not in protein_subunit_info_dic:
            continue

        #skip the case that gene name contain kda digit
        if ki in protein_kda_dic:
            if protein_kda_dic[ki][0]:
                gene_kda_pattern='\D'+protein_kda_dic[ki][1]+'(\D|$)'
                complie_gene_kda_pattern=re.compile(gene_kda_pattern)
                sr_gene_kda=complie_gene_kda_pattern.search(gene_name)
                if sr_gene_kda:
                    continue

        if ki in protein_subunit_info_dic:
            if len(protein_subunit_info_dic[ki][0])>0:
                complie_gene_greek=re.compile(gene_greek_pattern)
                sr_gene_greek=complie_gene_greek.search(gene_name)
                if sr_gene_greek:
                    #if the matched gene part is the same with the match protein part
                    #for example A <==> alpha
                    if sr_gene_greek.group(0)[0].upper()==protein_subunit_info_dic[ki][0][0].upper():
                        if ki in complex_dic:
                            complex_dic[ki].append([ki,gene_name,gene_name[:sr_gene_greek.span()[0]],'gene'])
                        else:
                            complex_dic[ki]=[[ki,gene_name,gene_name[:sr_gene_greek.span()[0]],'gene']]
                        #print('Only greek',[ki,gene_name,gene_name[:sr_gene_greek.span()[0]]])

            if len(protein_subunit_info_dic[ki][1])>0:
                complie_gene_greek=re.compile(gene_greek_pattern)
                sr_gene_greek=complie_gene_greek.search(gene_name)
                if sr_gene_greek:
                    #if the matched gene part is the same with the match protein part
                    #for example A1 <==> alpha-1
                    if sr_gene_greek.group(0).endswith(protein_subunit_info_dic[ki][1]) or \
                        sr_gene_greek.group(0).endswith(protein_subunit_info_dic[ki][1][0]):
                        if ki in complex_dic:
                            complex_dic[ki].append([ki,gene_name,gene_name[:sr_gene_greek.span()[0]],'gene'])
                        else:
                            complex_dic[ki]=[[ki,gene_name,gene_name[:sr_gene_greek.span()[0]],'gene']]
                        #print('Greek+digit:',[ki,gene_name,gene_name[:sr_gene_greek.span()[0]]])

            if len(protein_subunit_info_dic[ki][2])>0:
                #please notice that the re pattern match one more non digit letter here
                complie_single_digit=re.compile('\D'+protein_subunit_info_dic[ki][2]+'$')
                sr_single_digit=complie_single_digit.search(gene_name)
                
                if sr_single_digit:
                    if ki in complex_dic:
                        complex_dic[ki].append([ki,gene_name,gene_name[:sr_single_digit.span()[0]+1],'gene'])
                    else:
                        complex_dic[ki]=[[ki,gene_name,gene_name[:sr_single_digit.span()[0]+1],'gene']]
                    #print('Single digit',[ki,gene_name,gene_name[:sr_single_digit.span()[0]]])

            if len(protein_subunit_info_dic[ki][3])>0:
                complie_gene_single_letter=re.compile(protein_subunit_info_dic[ki][3]+'$')
                sr_gene_single_letter=complie_gene_single_letter.search(gene_name)
                #if ki=='Q7Z2R4':
                #    print(sr_gene_single_letter)
                if sr_gene_single_letter:
                    if ki in complex_dic:
                        complex_dic[ki].append([ki,gene_name,gene_name[:sr_gene_single_letter.span()[0]],'gene'])
                    else:
                        complex_dic[ki]=[[ki,gene_name,gene_name[:sr_gene_single_letter.span()[0]],'gene']]
                    #print('Single letter:',[ki,gene_name,gene_name[:sr_gene_single_letter.span()[0]]])

if uniprot_id in complex_dic:
    print(complex_dic[uniprot_id])

#here we can manually add some alternative name for a specific complex name 
#for example "Follicle-stimulating hormone" is a complex name, we add "Follicle stimulating hormone"
manual_add_complex_dic={'Follicle-stimulating hormone':'Follicle stimulating hormone',
                        }
for ki in complex_dic:
    for ei in complex_dic[ki]:
        if ei[2] in manual_add_complex_dic:
            manual_add_complex=[ei[0],manual_add_complex_dic[ei[2]]+ei[1][len(ei[2]):],manual_add_complex_dic[ei[2]]]
            complex_dic[ki].append(manual_add_complex)
#here also add the specse prefix for the gene name FSH-->hFSH
for ki in complex_dic:
    complex_name_of_this_id=[ei[2] for ei in complex_dic[ki]]
    for ei in complex_dic[ki]:
        if len(ei)==4:
            if not ei[2].startswith('h') and not ei[2].startswith('rh'):
                h_gene='h'+ei[2]
                rh_gene='rh'+ei[2]
                if h_gene not in complex_name_of_this_id:
                    complex_with_species_prefix=[ei[0],'h'+ei[1],h_gene,'gene']
                    complex_dic[ki].append(complex_with_species_prefix)
                if rh_gene not in complex_name_of_this_id:
                    complex_with_species_prefix=[ei[0],'rh'+ei[1],rh_gene,'gene']
                    complex_dic[ki].append(complex_with_species_prefix)

complex_file='./Glygen/OGER/uniprot_human_sprot_complex.csv'
with codecs.open(complex_file, 'w',encoding='utf-8') as tsvfile1:
    spamwriter = csv.writer(tsvfile1, delimiter='\t', quotechar='|')
    col_title=['cui','resource','original_id','term','preferred_term','entity_type']
    spamwriter.writerow(col_title)
    
    all_complex_list=[]
    for ki in complex_dic:
        for ei in complex_dic[ki]:
            row=['CUI-less','Swiss-Prot','','','','gene/protein']
            row[2]=ei[0]
            row[3]=ei[2]
            row[4]=ei[2]
            if row[3] not in all_complex_list:
                all_complex_list.append(row[3])
                spamwriter.writerow(row)
            


complex_file='./Glygen/OGER/complex_dic.csv'
with codecs.open(complex_file, 'w',encoding='utf-8') as tsvfile1:
    spamwriter = csv.writer(tsvfile1, delimiter='\t', quotechar='|')
    col_title=['complex','original name','type','uniprot_id']
    spamwriter.writerow(col_title)
    
    all_complex_list=[]
    for ki in complex_dic:
        #in complex_dic (uniprot id, original name, complex name)
        for ei in complex_dic[ki]:
            row=['']*len(col_title)
            row[0]=ei[2]
            row[1]=ei[1]
            row[2]=ei[1][len(ei[2]):]
            row[3]=ei[0]
            if row not in all_complex_list:
                all_complex_list.append(row)
                spamwriter.writerow(row)

#build a complex dictionary and dump it into json file
#then we can directly load it for future use
all_complex_list=[]
complex_dic_for_json={}
for ki in complex_dic:
    #in complex_dic (uniprot id, original name, complex name)
    for ei in complex_dic[ki]:
        if ei[2].lower() not in complex_dic_for_json:
            #in the dic, the key is complex name
            #the value is (type, uniprot id)
            complex_dic_for_json[ei[2].lower()]=[(ei[1][len(ei[2]):],ei[0])]
        else:
            if (ei[1][len(ei[2]):],ei[0]) not in complex_dic_for_json[ei[2].lower()]:
                complex_dic_for_json[ei[2].lower()].append((ei[1][len(ei[2]):],ei[0]))
                

complex_dic_json='./Glygen/OGER/complex_dic.json'
with open(complex_dic_json,'w') as dicfile:
        json.dump(complex_dic_for_json,dicfile)