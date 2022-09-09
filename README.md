# LLPSWise
computational prediction of biological liquid-liquid phase separation systems

https://www.biorxiv.org/content/10.1101/2022.07.25.501404v1

## Getting Started
1. pre-requirements:
see Code/packages.ini
2. data requirements:

a) Download PPI data （"BIOGRID-ORGANISM-Homo_sapiens-4.4.207.tab3.txt"） from biogrid

https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-4.4.211/

and change the path in main.py 

b) Create uniprot database from https://www.uniprot.org/help/downloads

``mkdir uniprotDB``

``wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz``

split the .dat file and rename those file as uniprotID.txt

3. start the query
python main.py --id [uniprotID]

4. check the results

change the path in app.py

``app = Flask(__name__, template_folder="../Code/condensateNetV2/templates")``

python app.py

then open a browser on the URL http://127.0.0.1:8000/llpswise/<targetid>
