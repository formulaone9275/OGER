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
                    protein_of_this_id=[pi[0] for pi in uniprot_kb_dic[row[2]]]
                    if row[3] not in protein_of_this_id:
                        spamwriter.writerow(row)
                    for li in uniprot_kb_dic[row[2]]:
                        row[3]=li[0]
                        row[-1]=li[1]
                        spamwriter.writerow(row)
                

def add_extra_items_in_dictionary(original_term_file):
    original_term_file='./Glygen/OGER/uniprot_human_sprot.csv'
    new_term_file='./Glygen/OGER/uniprot_human_sprot_added.csv'

    row_template=['CUI-less', 'Swiss-Prot', '', '', '', '']
    written_id=[]
    manual_list=[
        ['AMPA receptor','P19493','protein'],
        ['Tamm-Horsfall','P07911','protein'],
        ['IgM','P01871','gene'],
        ['IgE','P01854','gene'],
        ]
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
            #for the cases of coagulation factor ...
            if row[3].startswith('coagulation') or row[3].startswith('Coagulation'):
                if row[2] in coagulation_factor_dic:
                    coagulation_factor_dic[row[2]].append(row[3])
                else:
                    coagulation_factor_dic[row[2]]=[row[3]]
            #for the cases of genes, need to add the prefix of 'h' and 'rh' for each gene
            elif row[-1]=='gene' and len(row[3].split(' '))==1:
                if row[2] in gene_dic:
                    gene_dic[row[2]].append(row[3])
                else:
                    gene_dic[row[2]]=[row[3]]
            #for the cases of proeins with length of 2 like 'interleukin 5'
            #generalized form is: one word + digit/capitalized letter
            #need to add the hyphen connected and space connected and one word form
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
            elif row[-1]=='protein' and len(row[3].split(' '))==2:
                space_index=row[3].find(' ')
                second_part_of_protein=row[3][space_index+1:].strip()
                sr_capital_letters=complie_capital_letters.search(second_part_of_protein)
                if second_part_of_protein.isdigit() or sr_capital_letters:
                    print(row)
                    print(second_part_of_protein)
                    print(sr_capital_letters)
                    if row[2] in protein_with_two_words_dic:
                        protein_with_two_words_dic[row[2]].append(row[3])
                    else:
                        protein_with_two_words_dic[row[2]]=[row[3]]
    with open(new_term_file, 'w') as csvfile2:
            spamwriter = csv.writer(csvfile2, delimiter='\t', quotechar='|')
            spamwriter.writerow(col_title)

if __name__=='__main__':
    #build_term_list_for_human_sprot()
    #build_term_list_for_human_sprot_complex_or_protein()

    term_file='./Glygen/OGER/uniprot_human_sprot.csv'
    add_extra_items_in_dictionary(term_file)