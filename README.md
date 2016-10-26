## SSCO
Find the genes involved in start-stop codon overlaping in prokaryotic genome

## Strategy
This pipeline is built to find the genes involved in start-stop codon overlaping in prokaryotic genome.   
The idea is that for each start-stop codon overlaping motif in the genome, we check that if the downstream is a real gene.  
If the downstream is a real gene, we find if there is an annotated upstream gene. If we cannot find the upstream gene annotated, find the maximum-length upstream ORF (Open Reading Frame) in intergenic regions and annotate the sequence by blast.

## Preparation
You have to install Python, Biopython and Blast to run this pipline.  
You can download and install Python, Biopython and Blast from the link below:

Python https://www.python.org/ or sudo apt-get install python (Ubuntu)   
Biopython http://biopython.org/ or sudo apt-get install python-biopython (Ubuntu)   
Blast ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
#    

The nr database for blastx is needed.   
nr database for blastx  ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz    
Remmenber to use the command line below to make database before you use it

    makeblastdb -in <path> -dbtype prot -out <path>
#    

The genbank file is needed for analysis    
For exmaple, NC_000913.gbk is the genbank file for Escherichia coli str. K-12 substr. MG1655, complete genome   
You can download genbank file of specific specie from: ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Bacteria/  
Or you can download all genbank files of bacterial from :  ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Bacteria/all.gbk.tar.gz

## Getting Started 

Use NC_000913.gbk as the test file, an example     
Please ensure that the gkb file to be analyzed is under current directory or you may modify the path in scripts
Please ensure that all the result file is in the same directory with gbk file.
#

####command line #1

    python candidate_find.py

Candidate_find.py find and analyze all genbank files under the current directory.   
Example for result file: NC_002127_candidate.txt   
The result file contains information of all start-stop codon overlaping motif whose downstream is a real gene
#    

####command line #2

    python upstream.py

upstream.py is to find the upstream gene information. If we cannot find the upstream gene annotated, find the maximum-length upstream ORF (Open Reading Frame) in intergenic regions   
Example for result file: 
NC_000913_annotated_upstream_gene  (information of the downstream and upstream gene)
NC_000913_unknown_upstream_gene.fa (FASTA format file, used for blast)
#   

####command line #3

    python infor.py

infor.py is to find the function information of the downstream and upstream annotated gene   
Example for result file: NC_000913_annotated_upstream_gene.info
#   

####command line #4

    nohup blastx -query NC_000913_unknown_upstream_gene.fa -db <path of nr> -out NC_000913_unknown_upstream_gene.xml -outfmt 5 -evalue 0.001 -num-threads <Int> &

Ensure the output file is in XML format (format 5)   
Example for result file: NC_000913_unknown_upstream_gene.xml 
#   

####command line #5

    python blast_result.py

blast_result.py is used to analyze the output file of Blast.   
Example for result file: NC_002695_unknown_upstream_gene.info   
The result file contains function information of the downstream and upstream unknown gene   
You can change the select line of E-value in the script
  
#     
Contact information of the author: xiyupeng@iastate.edu
