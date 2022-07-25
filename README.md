# LLPSWise
computational prediction of biological liquid-liquid phase separation systems

## Getting Started
1. pre-requirements:
see Code/packages.ini
2. create uniprot database from https://www.uniprot.org/help/downloads
``mkdir uniprotDB``
``wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz``

split the .dat file and rename those file as uniprotID.txt

3. start the query
python main.py --id <uniprotID>

## Directory Structure
   deepchembed (master)
|--Data  
   |-- proteomeLLPSscore
|--Code  
   |-- app.py                       
   |-- main.py
   |-- model.py
   |-- packages.ini
   |-- uniprot_function3.py
