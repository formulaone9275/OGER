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

    with open(original_term_file, 'r') as csvfile1:
        spamreader = csv.reader(csvfile1, delimiter='\t', quotechar='|')
        with open(new_term_file, 'w') as csvfile2:
            spamwriter = csv.writer(csvfile2, delimiter='\t', quotechar='|')
            for row in spamreader:
                col_title=row
                spamwriter.writerow(col_title)
                break
            for row in spamreader:
                if row[2] in uniprot_kb_dic.keys():
                    if row[3] not in uniprot_kb_dic[row[2]]:
                        spamwriter.writerow(row)
                    for li in uniprot_kb_dic[row[2]]:
                        row[3]=li
                        spamwriter.writerow(row)
                

def build_term_list_for_human_sprot_complex():



if __name__=='__main__':
    build_term_list_for_human_sprot_complex()