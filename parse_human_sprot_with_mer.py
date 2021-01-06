from __future__ import print_function,division
import re,codecs,random,csv,json

uniprot_kb_file='./Glygen/OGER/uniprot_sprot_human.dat'

hmer_pattern='([Hh]\w*mer)\s?(;|of|with)'
com_hmer=re.compile(hmer_pattern)

uniprot_kb_dic={}
uniprot_kb_hmer_dic={}
with codecs.open(uniprot_kb_file,'r','utf8') as kbf:
    count=0
    for line in kbf:
        line=line.strip()

        if len(line)==0:
            continue
        if line.startswith('AC   '):
            id_list=line[3:].split(';')
            id_list=[i.strip() for i in id_list if i!='']
            for idi in id_list:
                uniprot_kb_dic[idi]=[]
            #print(id_list)
        elif line.startswith('DE   RecName:') or \
            line.startswith('DE   AltName: Full=') or \
                line.startswith('DE            Short='):
            #sometimes, there will be a section in {}
            parenthesis_left_ind=line.find('{')
            parenthesis_right_ind=line.find('}')
            if parenthesis_left_ind>0 and parenthesis_right_ind>0:
                line=line[:parenthesis_left_ind]+line[parenthesis_right_ind+1:]
            line_list=line.split('=')
            line_list=[i.strip() for i in line_list if i!='']
            #print(line_list)
            name_list=line_list[-1].split(';')
            name_list=[i.strip() for i in name_list if i!='']
            #print(name_list)
            for idi in id_list:
                uniprot_kb_dic[idi]+=name_list

        elif line.startswith('GN   Name=') or \
            line.startswith('GN   Synonyms='):
            #sometimes, there will be a section in {}
            parenthesis_left_ind=line.find('{')
            parenthesis_right_ind=line.find('}')
            if parenthesis_left_ind>0 and parenthesis_right_ind>0:
                line=line[:parenthesis_left_ind]+line[parenthesis_right_ind+1:]
            line_list=line[3:].split(';')
            line_list=[i.strip() for i in line_list if i!='']
            #print(line_list)
            name_list=[i.split('=')[-1] for i in line_list]
            name_list=[i.strip() for i in name_list]
            #in some case, there are some comma sparated gene name
            new_name_list=[]
            for ni in name_list:
                if ',' in ni:
                    gene_names=ni.split(',')
                    gene_names=[i.strip() for i in gene_names if i!='']
                    new_name_list+=gene_names
                else:
                    new_name_list.append(ni)
                    
            #name_list=[i.split(';')[0] for i in name_list]
            name_list=[i.strip() for i in new_name_list if i!='']
            #print(name_list)
            for idi in id_list:
                uniprot_kb_dic[idi]+=name_list
        elif line.startswith('CC   -!- SUBUNIT:'):
            sr=com_hmer.search(line[len('CC   -!- SUBUNIT:'):])
            if sr:
                print(line)
                for idi in id_list:
                    if idi in uniprot_kb_hmer_dic:
                        uniprot_kb_hmer_dic[idi]+=uniprot_kb_dic[idi]
                    else:
                        uniprot_kb_hmer_dic[idi]=uniprot_kb_dic[idi]
                    
json_file='./Glygen/OGER/uniprot_human_sprot_hmer.json'
with open(json_file,'w') as outfile:
    json.dump(uniprot_kb_hmer_dic,outfile)