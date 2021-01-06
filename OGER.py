from oger.ctrl.router import Router, PipelineServer
from nltk.corpus import words
import json,csv
conf = Router(termlist_path='./Glygen/OGER/uniprot_human_sprot.csv')
pl = PipelineServer(conf)

english_dict=words.words()
word_extra=['has']
for wi in word_extra:
    english_dict.append(wi)
amino_acid_shrot=['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Glu', 'Gln', 'Gly', 'His', 'Hyp', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Glp', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val',]
#pmid_list=[24342833,24036510,23285087,23242014,22160680,21712440,16679516,16040958,15504740]
#pmid_list=[9689040]
pmid_file='./Glygen/OGER/glygen_set.txt'
pmid_set=set()
with open(pmid_file, 'r') as f:
    for line in f:
        pmid_set.add(line.strip())
pmid_list=list(pmid_set)
pmid_str_list=[str(pi) for pi in pmid_list]
entity_dic={}

coll = pl.load_one(pmid_str_list, fmt='pubmed')
pl.process(coll)

protein_common_word_file='protein_common_word_file.csv'
with open(protein_common_word_file, 'w') as csvfile2:
    spamwriter = csv.writer(csvfile2, delimiter='\t', quotechar='|')
    spamwriter.writerow(['pmid','protein','protein id','start','end'])
    for pi in range(len(pmid_list)):
        print('pmid: ',pmid_list[pi])
        entities_list=[]
        for ei in coll[pi].iter_entities():
            if ei.info[2]=='Swiss-Prot' and ei.text not in english_dict and ei.text not in amino_acid_shrot:
                print((ei.text, ei.info[3], ei.start, ei.end))
                entities_list.append([ei.text,ei.info[3],int(ei.start),int(ei.end)])
                #print(ei.info)
            elif ei.info[2]=='Swiss-Prot' and ei.text in english_dict:
                spamwriter.writerow([str(pmid_list[pi]),ei.text,ei.info[3],int(ei.start),int(ei.end)])
        entity_dic[pmid_list[pi]]=entities_list

json_file='oger_entity.json'
with open(json_file,'w') as outfile:
    json.dump(entity_dic,outfile)
