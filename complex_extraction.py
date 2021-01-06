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
#for kDa match
kDa_pattern='\d+\s?kDa'

#for the information after subunit
after_subunit_alpha='[Aa]lpha'
after_subunit_beta='[B]eta'
after_subunit_single_digit='\d+'
after_subunit_single_letter='[A-Z]'

complie_subunit_pattern_1=re.compile(subunit_pattern_1)
complie_subunit_pattern_2=re.compile(subunit_pattern_2)

complie_kDa_pattern=re.compile(kDa_pattern)

complie_after_subunit_alpha=re.compile(after_subunit_alpha)
complie_after_subunit_beta=re.compile(after_subunit_beta)
complie_after_subunit_single_digit=re.compile(after_subunit_single_digit)
complie_after_subunit_single_letter=re.compile(after_subunit_single_letter)

protein_kda_dic={}
protein_subunit_info_dic={}
for ki in uniprot_kb_protein_dic.keys():
    for protein_name in uniprot_kb_protein_dic[ki]:

        sr_subunit_1=complie_subunit_pattern_1.search(protein_name)
        sr_subunit_2=complie_subunit_pattern_2.search(protein_name)
        sr_kda=complie_kDa_pattern.search(protein_name)

        #find the kda digit to remove relevant the gene name
        #[True/False, kDa digit]
        protein_kda_dic[ki]=[False,'']
        #[alpha, beta, digit, single letter]
        protein_subunit_info_dic[ki]=['','','','']
        if sr_kda:

            kda_digit=protein_name[sr_kda.span()[0]:sr_kda.span()[1]-3]
            kda_digit=kda_digit.strip()
            protein_kda_dic[ki]=[True,kda_digit]

        if sr_subunit_1:
            if ki in complex_dic:
                complex_dic[ki].append([ki,protein_name,protein_name[:sr_subunit_1.span()[0]]])
            else:
                complex_dic[ki]=[[ki,protein_name,protein_name[:sr_subunit_1.span()[0]]]]
            #print('ID:',ki)
            #print('Original protein name: ',protein_name)
            #print(sr_subunit_1.span())
            #print('Complex name: ',protein_name[:sr_subunit_1.span()[0]])
            after_subunit_str=protein_name[sr_subunit_1.span()[1]+1:]
            after_subunit_str=after_subunit_str.strip()
            sr_alpha=complie_after_subunit_alpha.search(after_subunit_str)
            sr_beta=complie_after_subunit_beta.search(after_subunit_str)
            if sr_alpha:
                protein_subunit_info_dic[ki][0]='alpha'
            if sr_beta:
                protein_subunit_info_dic[ki][1]='beta'
            if after_subunit_str.isdigit():
                protein_subunit_info_dic[ki][2]=after_subunit_str
            if len(after_subunit_str)==1:
                protein_subunit_info_dic[ki][3]=after_subunit_str

        elif sr_subunit_2:
            if ki in complex_dic:
                complex_dic[ki].append([ki,protein_name,protein_name[:sr_subunit_2.span()[0]]])
            else:
                complex_dic[ki]=[[ki,protein_name,protein_name[:sr_subunit_2.span()[0]]]]
            #print('ID:',ki)
            #print('Original protein name: ',protein_name)
            #print(sr_subunit_2.span())
            #print('Complex name: ',protein_name[:sr_subunit_2.span()[0]])
            after_subunit_str=protein_name[sr_subunit_2.span()[1]+1:]
            after_subunit_str=after_subunit_str.strip()
            sr_alpha=complie_after_subunit_alpha.search(after_subunit_str)
            sr_beta=complie_after_subunit_beta.search(after_subunit_str)
            if sr_alpha:
                protein_subunit_info_dic[ki][0]='alpha'
            if sr_beta:
                protein_subunit_info_dic[ki][1]='beta'
            if after_subunit_str.isdigit():
                protein_subunit_info_dic[ki][2]=after_subunit_str
            if len(after_subunit_str)==1:
                protein_subunit_info_dic[ki][3]=after_subunit_str


#for the information after subunit
gene_alpha_pattern='A\d*$'
gene_beta_pattern='B\d*$'
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
                complie_gene_alpha=re.compile(gene_alpha_pattern)
                sr_gene_alpha=complie_gene_alpha.search(gene_name)
                if sr_gene_alpha:
                    if ki in complex_dic:
                        complex_dic[ki].append([ki,gene_name,gene_name[:sr_gene_alpha.span()[0]]])
                    else:
                        complex_dic[ki]=[[ki,gene_name,gene_name[:sr_gene_alpha.span()[0]]]]
                    print([ki,gene_name,gene_name[:sr_gene_alpha.span()[0]]])

            elif len(protein_subunit_info_dic[ki][1])>0:
                complie_gene_beta=re.compile(gene_beta_pattern)
                sr_gene_beta=complie_gene_beta.search(gene_name)
                if sr_gene_beta:
                    if ki in complex_dic:
                        complex_dic[ki].append([ki,gene_name,gene_name[:sr_gene_beta.span()[0]]])
                    else:
                        complex_dic[ki]=[[ki,gene_name,gene_name[:sr_gene_beta.span()[0]]]]
                    print([ki,gene_name,gene_name[:sr_gene_beta.span()[0]]])

            elif len(protein_subunit_info_dic[ki][2])>0:
                complie_single_digit=re.compile(protein_subunit_info_dic[ki][2]+'$')
                sr_single_digit=complie_single_digit.search(gene_name)
                if sr_single_digit:
                    if ki in complex_dic:
                        complex_dic[ki].append([ki,gene_name,gene_name[:sr_single_digit.span()[0]]])
                    else:
                        complex_dic[ki]=[[ki,gene_name,gene_name[:sr_single_digit.span()[0]]]]
                    print([ki,gene_name,gene_name[:sr_single_digit.span()[0]]])

            elif len(protein_subunit_info_dic[ki][3])>0:
                complie_gene_single_letter=re.compile(protein_subunit_info_dic[ki][3]+'$')
                sr_gene_single_letter=complie_gene_single_letter.search(gene_name)
                if sr_gene_single_letter:
                    if ki in complex_dic:
                        complex_dic[ki].append([ki,gene_name,gene_name[:sr_gene_single_letter.span()[0]]])
                    else:
                        complex_dic[ki]=[[ki,gene_name,gene_name[:sr_gene_single_letter.span()[0]]]]
                    print([ki,gene_name,gene_name[:sr_gene_single_letter.span()[0]]])

complex_file='./Glygen/OGER/complex_dic.csv'
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
            if row not in all_complex_list:
                all_complex_list.append(row)
                spamwriter.writerow(row)