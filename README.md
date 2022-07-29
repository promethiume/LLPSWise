# LLPSWise
computational prediction of biological liquid-liquid phase separation systems

## Getting Started
1. pre-requirements:
see Code/packages.ini
2. data requirements:

a) Download PPI data from biogrid

https://downloads.thebiogrid.org/BioGRID

Or use our version(already uploaded in this repository):

"BIOGRID-ORGANISM-Homo_sapiens-4.4.207.tab3.txt"

b) Create uniprot database from https://www.uniprot.org/help/downloads

``mkdir uniprotDB``

``wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz``

split the .dat file and rename those file as uniprotID.txt

3. start the query
python main.py --id [uniprotID]
