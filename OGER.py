from oger.ctrl.router import Router, PipelineServer
from nltk.corpus import words
import json,csv,re


def protein_detection_OGER(term_file,output_json,output_csv):
    conf = Router(termlist_path=term_file)
    pl = PipelineServer(conf)

    english_dict=words.words()
    word_extra=['has','Post','To']
    for wi in word_extra:
        english_dict.append(wi)
    amino_acid_shrot=['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Glu', 'Gln', 'Gly', 'His', 'Hyp', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Glp', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val',]
    #pmid_list=[24342833,24036510,23285087,23242014,22160680,21712440,16679516,16040958,15504740]
    #pmid_list=[17286803]
    
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
    #print(coll[0][0].text)
    filter_out_protein_pattern='^([A-Za-z]\s?\d|[A-Za-z]{2})$'
    com_protein_filter_out=re.compile(filter_out_protein_pattern)
    filter_out_protein_list=['and 1','at 1','solved at 1','GP1','GP2']
    protein_common_word_file='protein_common_word_file.csv'
    with open(output_csv, 'w') as csvfile1:
        spamwriter1 = csv.writer(csvfile1, delimiter='\t', quotechar='|')
        spamwriter1.writerow(['pmid','protein/complex','protein id','start','end'])
        with open(protein_common_word_file, 'w') as csvfile2:
            spamwriter2 = csv.writer(csvfile2, delimiter='\t', quotechar='|')
            spamwriter2.writerow(['pmid','protein','protein id','start','end'])
            for pi in range(len(pmid_list)):
                #print('pmid: ',pmid_list[pi])
                entities_list=[]
                for ei in coll[pi].iter_entities():
                    #print(ei.text)
                    sr=com_protein_filter_out.search(ei.text)
                    name_is_digit=ei.text.isdigit()
                    
                    if ei.info[2]=='Swiss-Prot' and ei.text not in english_dict and \
                         ei.text not in amino_acid_shrot and not sr and not name_is_digit and \
                             ei.text not in filter_out_protein_list:
                        #print((ei.text, ei.info[3], ei.start, ei.end))
                        entities_list.append([ei.text,ei.info[3],int(ei.start),int(ei.end)])
                        spamwriter1.writerow([str(pmid_list[pi]),ei.text,ei.info[3],int(ei.start),int(ei.end)])
                        #print(ei.info)
                    elif ei.info[2]=='Swiss-Prot' and ei.text in english_dict:
                        spamwriter2.writerow([str(pmid_list[pi]),ei.text,ei.info[3],int(ei.start),int(ei.end)])
                #in some cases, OGER will detect entities that are overlapped to each other
                #like ['Tissue plasminogen activator', 'P00750', 0, 28], ['plasminogen', 'P00747', 7, 18]
                
                entities_list_non_overlap=[]
                for entity_i in entities_list:
                    overlap_sign=False
                    for entity_ii in entities_list:
                        if (entity_i[2]>entity_ii[2] and entity_i[3]<=entity_ii[3]) or \
                        (entity_i[2]>=entity_ii[2] and entity_i[3]<entity_ii[3]) or \
                            (entity_i[2]>entity_ii[2] and entity_i[3]<entity_ii[3]):
                            overlap_sign=True
                            break
                    if not overlap_sign:
                        entities_list_non_overlap.append(entity_i)
                entity_dic[pmid_list[pi]]=entities_list_non_overlap
    #print(entity_dic['12654314'])
    with open(output_json,'w') as outfile:
        json.dump(entity_dic,outfile)


if __name__=='__main__':

    
    
    
    term_file='./Glygen/OGER/uniprot_human_sprot_final.csv'
    output_json='./Glygen/OGER/oger_entity_sorted.json'
    output_csv='./Glygen/OGER/oger_entity_sorted.csv'
    protein_detection_OGER(term_file,output_json,output_csv)
    
    '''
    term_file='./Glygen/OGER/uniprot_human_sprot_cp_sorted.csv'
    output_json='./Glygen/OGER/oger_entity_cp_sorted.json'
    output_csv='./Glygen/OGER/oger_entity_cp_sorted.csv'
    protein_detection_OGER(term_file,output_json,output_csv)

    term_file='./Glygen/OGER/uniprot_human_sprot_cp.csv'
    output_json='./Glygen/OGER/oger_entity_cp.json'
    output_csv='./Glygen/OGER/oger_entity_cp.csv'
    protein_detection_OGER(term_file,output_json,output_csv)


    term_file='./Glygen/OGER/uniprot_human_sprot_complex_sorted.csv'
    output_json='./Glygen/OGER/oger_entity_complex_sorted.json'
    output_csv='./Glygen/OGER/oger_entity_complex_sorted.csv'
    protein_detection_OGER(term_file,output_json,output_csv)

    '''
    
    