from __future__ import print_function,division
import re,codecs,random,csv,json
import numpy as np
import pandas as pd


def build_term_list_for_human_sprot():
    original_term_file='./Glygen/OGER/uniprot.csv'
    new_term_file='./Glygen/OGER/uniprot_human_sprot.csv'

    json_file='./Glygen/OGER/uniprot_human_sprot.json'
    with open(json_file) as jfile:
        uniprot_kb_dic=json.load(jfile)

    written_id=[]
    with open(original_term_file, 'r') as csvfile1:
        spamreader = csv.reader(csvfile1, delimiter='\t', quotechar='|')
        with open(new_term_file, 'w') as csvfile2:
            spamwriter = csv.writer(csvfile2, delimiter='\t', quotechar='|')
            for row in spamreader:
                col_title=row
                spamwriter.writerow(col_title)
                break
            for row in spamreader:
                if row[2] in uniprot_kb_dic.keys() and row[2] not in written_id:
                    written_id.append(row[2])
                    protein_of_this_id=[pi[0] for pi in uniprot_kb_dic[row[2]]]
                    if row[3] not in protein_of_this_id:
                        spamwriter.writerow(row)
                    for li in uniprot_kb_dic[row[2]]:
                        row[3]=li[0]
                        row[-1]=li[1]
                        spamwriter.writerow(row)
                

def build_term_list_for_human_sprot_complex():
    print('Already implemented it in the complex_extraction.py!')


def build_term_list_for_human_sprot_complex_or_protein():
    original_term_file='./Glygen/OGER/uniprot.csv'
    new_term_file='./Glygen/OGER/uniprot_human_sprot_cp.csv'

    json_file='./Glygen/OGER/uniprot_human_sprot_hmer.json'
    with open(json_file) as jfile:
        uniprot_kb_dic=json.load(jfile)
    written_id=[]
    with open(original_term_file, 'r') as csvfile1:
        spamreader = csv.reader(csvfile1, delimiter='\t', quotechar='|')
        with open(new_term_file, 'w') as csvfile2:
            spamwriter = csv.writer(csvfile2, delimiter='\t', quotechar='|')
            for row in spamreader:
                col_title=row
                spamwriter.writerow(col_title)
                break
            for row in spamreader:
                if row[2] in uniprot_kb_dic.keys() and row[2] not in written_id:
                    written_id.append(row[2])
                    protein_of_this_id=[pi for pi in uniprot_kb_dic[row[2]]]
                    if row[3] not in protein_of_this_id:
                        spamwriter.writerow(row)
                    for li in uniprot_kb_dic[row[2]]:
                        row[3]=li
                        row[-2]=li
                        spamwriter.writerow(row)
                
def filter_out_protein_in_dictionary(original_term_file,new_term_file):

    protein_pattern='^([A-Za-z]\s?\d|[A-Za-z]{2})$'
    com_protein=re.compile(protein_pattern)

    delete_protein_list=[('GP130','P42704'),('GP1','O00178'),('GP2','P55259')]
    with open(new_term_file, 'w') as csvfile2:
            spamwriter = csv.writer(csvfile2, delimiter='\t', quotechar='|')
            with open(original_term_file, 'r') as csvfile1:
                spamreader = csv.reader(csvfile1, delimiter='\t', quotechar='|')
                
                for row in spamreader:
                    col_title=row
                    break
                spamwriter.writerow(col_title)
                for row in spamreader:
                    sr=com_protein.search(row[3])
                    if sr:
                        continue
                    need_to_delete_sign=False
                    for dpi in delete_protein_list:
                        if row[3]==dpi[0] and row[2]==dpi[1]:
                            need_to_delete_sign=True
                            break
                    if need_to_delete_sign:
                        continue
                    spamwriter.writerow(row)

def add_extra_items_in_dictionary(original_term_file):
    original_term_file='./Glygen/OGER/uniprot_human_sprot.csv'
    new_term_file='./Glygen/OGER/uniprot_human_sprot_added.csv'
    json_file='./Glygen/OGER/uniprot_human_sprot.json'
    with open(json_file) as jfile:
        uniprot_kb_dic=json.load(jfile)

    row_template=('CUI-less', 'Swiss-Prot', '', '', '', '')
    protein_gene_extra_list=[]
    original_gene_protein_list=[]
    manual_list=[]
    with open('./Glygen/OGER/manual_add_protein.txt', 'r') as mfile:
        for line in mfile:
            line=line.strip()
            line_split=line.split('\t')
            line_split=[li for li in line_split if li!='']
            #some ids are given in lower case
            line_split[-1]=line_split[-1].upper()
            line_split.append('protein')
            manual_list.append(line_split)
    print('Manually added list: ',manual_list)
    coagulation_factor_dic={}
    gene_dic={}
    protein_with_hyphen_dic={}
    protein_with_two_words_dic={}
    capital_letters_pattern='^[A-Z]+$'
    complie_capital_letters=re.compile(capital_letters_pattern)
    with open(original_term_file, 'r') as csvfile1:
        spamreader = csv.reader(csvfile1, delimiter='\t', quotechar='|')
        
        for row in spamreader:
            col_title=row
            break
        for row in spamreader:
            original_gene_protein_list.append(row)
            #case 1: for the cases of coagulation factor ...
            if row[3].startswith('coagulation') or row[3].startswith('Coagulation'):
                if row[2] in coagulation_factor_dic:
                    coagulation_factor_dic[row[2]].append(row[3])
                else:
                    coagulation_factor_dic[row[2]]=[row[3]]
            #case 2: for the cases of genes, need to add the prefix of 'h' and 'rh' for each gene
            elif row[-1]=='gene' and len(row[3].split(' '))==1:
                if row[2] in gene_dic:
                    gene_dic[row[2]].append(row[3])
                else:
                    gene_dic[row[2]]=[row[3]]
            #case 3: for the cases of proeins with hyphen in them like 'interleukin-5'
            #need to add the immediately connected one word form and space connected and two words form
            elif row[-1]=='protein' and len(row[3].split(' '))==1 and '-' in row[3]:
                hyphen_index=row[3].find('-')
                second_part_of_protein=row[3][hyphen_index+1:].strip()
                sr_capital_letters=complie_capital_letters.search(second_part_of_protein)
                if second_part_of_protein.isdigit() or sr_capital_letters:
                    #print(row)
                    #print(second_part_of_protein)
                    #print(sr_capital_letters)
                    if row[2] in protein_with_hyphen_dic:
                        protein_with_hyphen_dic[row[2]].append(row[3])
                    else:
                        protein_with_hyphen_dic[row[2]]=[row[3]]
            #case 4: for the cases of proeins with length of 2 like 'interleukin 5'
            #generalized form is: one word + digit/capitalized letter
            #need to add the hyphen connected one word form and immediately connected one word form
            elif row[-1]=='protein' and len(row[3].split(' '))==2:
                space_index=row[3].find(' ')
                second_part_of_protein=row[3][space_index+1:].strip()
                sr_capital_letters=complie_capital_letters.search(second_part_of_protein)
                if second_part_of_protein.isdigit() or sr_capital_letters:
                    #print(row)
                    #print(second_part_of_protein)
                    #print(sr_capital_letters)
                    if row[2] in protein_with_two_words_dic:
                        protein_with_two_words_dic[row[2]].append(row[3])
                    else:
                        protein_with_two_words_dic[row[2]]=[row[3]]
    #generate the rows for the extra protein names
    #first the the items in the manual list
    for mii in manual_list:
        row_add=list(row_template)
        row_add[2]=mii[1]
        row_add[3]=mii[0]
        row_add[4]=mii[0]
        row_add[5]=mii[-1]
        protein_gene_extra_list.append(row_add)
        print(row_add)
    #for case 1
    for ki in coagulation_factor_dic.keys():
        protein_of_this_id=[pi[0] for pi in uniprot_kb_dic[ki]]
        for pi in coagulation_factor_dic[ki]:
            protein_wo_coagulation=pi[len('coagulation '):]
            protein_wo_coagulation_split=protein_wo_coagulation.split(' ')
            
            if protein_wo_coagulation not in protein_of_this_id:
                row_add=list(row_template)
                row_add[2]=ki
                row_add[3]=protein_wo_coagulation
                row_add[4]=pi
                row_add[5]='protein'
                protein_gene_extra_list.append(row_add)
                #print(row_add)
            if len(protein_wo_coagulation_split)==2:
                row_add=list(row_template)
                row_add[2]=ki
                row_add[3]='f'+protein_wo_coagulation_split[1]
                row_add[4]=pi
                row_add[5]='protein'
                protein_gene_extra_list.append(row_add)
                #print(row_add)
    #for case 2
    for ki in gene_dic.keys():
        protein_of_this_id=[pi[0] for pi in uniprot_kb_dic[ki]]
        for pi in gene_dic[ki]:
            if not pi.startswith('h') and not pi.startswith('rh'):
                h_gene='h'+pi
                if h_gene not in protein_of_this_id:
                    row_add=list(row_template)
                    row_add[2]=ki
                    row_add[3]=h_gene
                    row_add[4]=pi
                    row_add[5]='gene'
                    protein_gene_extra_list.append(row_add)
                    #print(row_add)
                rh_gene='rh'+pi
                if rh_gene not in protein_of_this_id:
                    row_add=list(row_template)
                    row_add[2]=ki
                    row_add[3]=rh_gene
                    row_add[4]=pi
                    row_add[5]='gene'
                    protein_gene_extra_list.append(row_add)
                    #print(row_add)
    #for case 3
    for ki in protein_with_hyphen_dic.keys():
        protein_of_this_id=[pi[0] for pi in uniprot_kb_dic[ki]]
        for pi in protein_with_hyphen_dic[ki]:
            protein_with_space=pi.replace('-',' ')
            protein_one_word=pi.replace('-','')
        
            if protein_with_space not in protein_of_this_id:
                row_add=list(row_template)
                row_add[2]=ki
                row_add[3]=protein_with_space
                row_add[4]=pi
                row_add[5]='protein'
                if row_add[3] not in ['And 1','AT 1']:
                    protein_gene_extra_list.append(row_add)
                #print(row_add)
            if protein_one_word not in protein_of_this_id:
                row_add=list(row_template)
                row_add[2]=ki
                row_add[3]=protein_one_word
                row_add[4]=pi
                row_add[5]='protein'
                protein_gene_extra_list.append(row_add)
                #print(row_add)
    #for case 4
    for ki in protein_with_two_words_dic.keys():
        protein_of_this_id=[pi[0] for pi in uniprot_kb_dic[ki]]
        for pi in protein_with_two_words_dic[ki]:
            protein_with_hyphen=pi.replace(' ','-')
            protein_one_word=pi.replace(' ','')
        
            if protein_with_hyphen not in protein_of_this_id:
                row_add=list(row_template)
                row_add[2]=ki
                row_add[3]=protein_with_hyphen
                row_add[4]=pi
                row_add[5]='protein'
                protein_gene_extra_list.append(row_add)
                #print(row_add)
            if protein_one_word not in protein_of_this_id:
                row_add=list(row_template)
                row_add[2]=ki
                row_add[3]=protein_one_word
                row_add[4]=pi
                row_add[5]='protein'
                protein_gene_extra_list.append(row_add)
                #print(row_add)
    with open(new_term_file, 'w') as csvfile2:
            spamwriter = csv.writer(csvfile2, delimiter='\t', quotechar='|')
            spamwriter.writerow(col_title)
            for ri in protein_gene_extra_list+original_gene_protein_list:
                spamwriter.writerow(ri)

def sort_term_file(term_file):
    term_dic={}
    with open(term_file, 'r') as csvfile1:
        spamreader = csv.reader(csvfile1, delimiter='\t', quotechar='|')
        
        for row in spamreader:
            col_title=row
            break
        for row in spamreader:
            if len(row[3])==0:
                continue
            term_dic[row[3]]=row
    sorted_term_dic={}
    for k in sorted(term_dic, key=len):
        sorted_term_dic[k] = term_dic[k]
    #print(sorted_term_dic.keys()[:10])
    #print(sorted_term_dic['Chorionic gonadotropin'])
    new_term_file='./Glygen/OGER/uniprot_human_sprot_sorted.csv'
    with open(new_term_file, 'w') as csvfile2:
            spamwriter = csv.writer(csvfile2, delimiter='\t', quotechar='|')
            spamwriter.writerow(col_title)
            for ri in sorted_term_dic:
                spamwriter.writerow(sorted_term_dic[ri])

if __name__=='__main__':
    #build_term_list_for_human_sprot()
    #build_term_list_for_human_sprot_complex_or_protein()

    term_file='./Glygen/OGER/uniprot_human_sprot.csv'
    add_extra_items_in_dictionary(term_file)
    term_file='./Glygen/OGER/uniprot_human_sprot_added.csv'
    sort_term_file(term_file)