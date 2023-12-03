# OGER for GlyGen
1. Raw file

The raw file (uniprot_sprot_human.dat) is either from Uniprot website or from Hongzhan (DBI);

2. Parse the raw file

The code to parse the raw file could be found at:parse_human_sprot.py;  

3. Build the dictionary

To build the dictionary for OGER, we could employ the code in build_term_list.py. As you can see from the code file, we build two dictionary files, one is regular (uniprot_human_sprot_final.csv), and other one is sorted (uniprot_human_sprot_sorted.csv). This is because that OGER will return different results using different dictionaries (the items in the dictionaries are the same, but in different order), and this is probably a bug of OGER. So we build two dictionaries and merge the OGER outputs.

4. Deployment

After we have the OGER dictionaries, we could deploy them in Docker using the code file OGER.py
