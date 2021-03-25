from oger.ctrl.router import Router, PipelineServer
from nltk.corpus import words
import json,csv,re,os


def protein_detection_OGER(term_file,full_length_txt_path):
    conf = Router(termlist_path=term_file)
    pl = PipelineServer(conf)

    english_dict=words.words()
    word_extra=['has','Post','To','Large']
    word_delete_from_dic=['thyroglobulin','albumin','insulin','rhodopsin',\
        'thyroglobulin','cholinesterase','prothrombin','prolactin','catalase','chymase','tyrosinase',\
            'elastin','renin','agrin',]
    for wi in word_extra:
        english_dict.append(wi)
    for wi in word_delete_from_dic:
        if wi in english_dict:
            english_dict.remove(wi)
    amino_acid_shrot=['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Glu', 'Gln', 'Gly', 'His', 'Hyp', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Glp', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val',]
    

    entity_dic={}
    #process txt file one by one
    for fi in os.listdir(full_length_txt_path):
        if not fi.endswith('txt'):
            continue
        
        coll = pl.load_one(full_length_txt_path+fi, 'txt')

        file_id=fi.split('.')[0]
        pl.process(coll)
        #print(coll[0][0].text)
        filter_out_protein_pattern='^([A-Za-z]\s?\d|[A-Za-z]{2})$'
        com_protein_filter_out=re.compile(filter_out_protein_pattern)
        filter_out_protein_list=['and 1','at 1','solved at 1','GP1','GP2','at 2','Fig']
    
        entities_list=[]
        
        for ei in coll[0].iter_entities():
            #print(ei.text)
            
            sr=com_protein_filter_out.search(ei.text)
            name_is_digit=ei.text.isdigit()
            
            if ei.info[2]=='Swiss-Prot' and ei.text not in english_dict and \
                    ei.text not in amino_acid_shrot and not sr and not name_is_digit and \
                        ei.text not in filter_out_protein_list:
                #print((ei.text, ei.info[3], ei.start, ei.end))
                entities_list.append([ei.text,ei.info[3],int(ei.start),int(ei.end)])
                
    
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
        entity_dic[file_id]=entities_list_non_overlap
            
    return entity_dic


def merge_oger_detections(term_file,term_file_sorted,output_json):

    #use different term file to detect entities and then merge them
    full_length_txt_path='./Glygen/OGER/full_length_txt/'
    entity_dic=protein_detection_OGER(term_file,full_length_txt_path)
    
    entity_dic_sorted=protein_detection_OGER(term_file_sorted,full_length_txt_path)

    #this is for logging the protein names that are normalized to more than one uniprot ids
    multiple_ids_protein_dic={}
    #for example in 11067851, use different term_file will give different entity detection results
    entity_dic_merged={}
    for ki in entity_dic.keys():
        entity_dic_merged[ki]=[]
        for vi in entity_dic[ki]:
            if vi not in entity_dic_merged[ki]:
                entity_dic_merged[ki].append(vi)
        for vi in entity_dic_sorted[ki]:
            if vi not in entity_dic_merged[ki]:
                entity_dic_merged[ki].append(vi)
        
        #for some entities like TSP-1['TSP-1', 'Q6UWB4', 90, 95], ['TSP-1', 'P07996', 90, 95], ['TSP-1', 'Q9BXG8', 90, 95]
        #different ids will be assigned to it, we need delete the duplicates
        #if the id is also detected in the same abstract, we will use that one
        #otherwise we will choose the first one
        deduplicate_entity_list=[]
        duplicated_entity_name_and_offset={}
        entity_name_and_offset=[]
        for ei in entity_dic_merged[ki]:
            if (ei[0],ei[2],ei[3]) not in entity_name_and_offset:
                entity_name_and_offset.append((ei[0],ei[2],ei[3]))
            else:
                duplicated_entity_name_and_offset[(ei[0],ei[2],ei[3])]=set()
                for eii in entity_dic_merged[ki]:
                    if (ei[0],ei[2],ei[3])==(eii[0],eii[2],eii[3]):
                        duplicated_entity_name_and_offset[(ei[0],ei[2],ei[3])].add(eii[1])
                    
        
        multiple_ids_protein_dic[ki]=duplicated_entity_name_and_offset

        for ei in entity_dic_merged[ki]:
            if (ei[0],ei[2],ei[3]) not in duplicated_entity_name_and_offset:
                deduplicate_entity_list.append(ei)
        for di in duplicated_entity_name_and_offset.keys():
            found_associated_id=False
            for ei in entity_dic_merged[ki]:
                if ei[0]!=di[0] and ei[1] in duplicated_entity_name_and_offset[di]:
                    found_associated_id=True
                    deduplicate_entity_list.append([di[0],ei[1],di[1],di[2]])
                    break
            if not found_associated_id:
                deduplicate_entity_list.append([di[0],list(duplicated_entity_name_and_offset[di])[0],di[1],di[2]])
        #["organic anion transporting polypeptide 1B1", "Q9Y6L6", 46, 88]
        #["organic anion transporting polypeptide 1", "P46721", 46, 86]
        entities_list_non_overlap=[]
        for entity_i in deduplicate_entity_list:
            overlap_sign=False
            for entity_ii in deduplicate_entity_list:
                if (entity_i[2]>entity_ii[2] and entity_i[3]<=entity_ii[3]) or \
                (entity_i[2]>=entity_ii[2] and entity_i[3]<entity_ii[3]) or \
                    (entity_i[2]>entity_ii[2] and entity_i[3]<entity_ii[3]):
                    overlap_sign=True
                    break
            if not overlap_sign:
                entities_list_non_overlap.append(entity_i)

        entity_dic_merged[ki]=entities_list_non_overlap

    #print(entity_dic_merged)
    with open(output_json,'w') as outfile:
        json.dump(entity_dic_merged,outfile)
    #this is for logging the protein names that are normalized to more than one uniprot ids
    row_list=[]
    for ki in multiple_ids_protein_dic.keys():
        for pi in multiple_ids_protein_dic[ki].keys():
            mpi_row=[ki,pi[0],','.join(list(multiple_ids_protein_dic[ki][pi]))]
            row_list.append(tuple(mpi_row))
    row_list=list(set(row_list))
    with open('proteins_with_multiple_normalized_ids.csv', 'w') as csvfile2:
            spamwriter_mpi = csv.writer(csvfile2, delimiter='\t', quotechar='|')
            spamwriter_mpi.writerow(['pmid','protein','protein ids'])
            for ri in row_list:
                spamwriter_mpi.writerow(list(ri))

if __name__=='__main__':
    '''
    term_file='./Glygen/OGER/uniprot_human_sprot_complex_final.csv'
    term_file_sorted='./Glygen/OGER/uniprot_human_sprot_complex_sorted.csv'
    output_json='./Glygen/OGER/glygen_large_oger_entity_complex.json'

    merge_oger_detections(term_file,term_file_sorted,output_json)
    
    '''
    #for proteins
    term_file='./Glygen/OGER/uniprot_human_sprot_final.csv'
    term_file_sorted='./Glygen/OGER/uniprot_human_sprot_sorted.csv'
    output_json='./Glygen/OGER/glygen_full_length_oger_entity.json'

    merge_oger_detections(term_file,term_file_sorted,output_json)
   
    '''
    #for complex-protein
    term_file='./Glygen/OGER/uniprot_human_sprot_cp_final.csv'
    term_file_sorted='./Glygen/OGER/uniprot_human_sprot_cp_sorted.csv'
    output_json='./Glygen/OGER/glygen_large_oger_entity_cp.json'

    merge_oger_detections(term_file,term_file_sorted,output_json)
    '''

    
    
    