from __future__ import print_function,division
import re,codecs,random,csv,json

uniprot_kb_file='./Glygen/OGER/uniprot_sprot.dat'

uniprot_kb_id_to_ac_dic={}
uniprot_kb_ac_to_id_dic={}
with codecs.open(uniprot_kb_file,'r','utf8') as kbf:
    count=0
    previous_line_is_ac=False
    for line in kbf:
        line=line.strip()

        if len(line)==0:
            continue
        if line.startswith('ID   '):
            previous_line_is_ac=False
            id_list=line[3:].split(' ')
            id_list=[i.strip() for i in id_list if i!='']
            id=id_list[0]
        elif line.startswith('AC   '):
            #if there are many AC's for one protein
            #only extract the first one and ignore the others
            if previous_line_is_ac:
                continue
            previous_line_is_ac=True
            ac_list=line[3:].split(';')
            ac_list=[i.strip() for i in ac_list if i!='']
            ac_list=ac_list[:1]
            uniprot_kb_id_to_ac_dic[id]=ac_list[0]
            uniprot_kb_ac_to_id_dic[ac_list[0]]=id
            print(id,'==>',ac_list[0])
        else:
            previous_line_is_ac=False

           

json_file='./Glygen/OGER/uniprot_kb_id_to_ac.json'
with open(json_file,'w') as outfile:
    json.dump(uniprot_kb_id_to_ac_dic,outfile)

json_file='./Glygen/OGER/uniprot_kb_ac_to_id.json'
with open(json_file,'w') as outfile:
    json.dump(uniprot_kb_ac_to_id_dic,outfile)