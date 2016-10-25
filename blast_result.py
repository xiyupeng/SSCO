#!usr/bin/python
#filename: blast_result.py
from Bio.Blast import NCBIXML
from Bio import SeqIO
import sys
import os
#read the result xml file from blast

path=os.getcwd()
E_VALUE=10**(-10)   # you can change the E-value here

origin=sys.stdout
for x in os.listdir(path):
	if ("_unknown_upstream_gene" in x and "xml" in x):
		specie=x.split('_unknown')[0]
		infor_file=str(x.split('.')[0]+".info")
		file=open(infor_file,'w')
		sys.stdout=file
		blast_file=open(x,'r')
		gb_file=str(specie+".gbk")
		for y in os.listdir(path):
			if gb_file in y:
				genename={}
				genefunction={}
				geneproduct={}
				gbfile=open(gb_file,"r")
				for gb_record in SeqIO.parse(gbfile, "genbank") :
					number=len(gb_record.features)
					for j in range(0,number):
						gb_feature = gb_record.features[j]
						genename[j]=str(gb_feature.qualifiers.get('gene')).replace ("['","").replace("']","")
						genefunction[j]=str(gb_feature.qualifiers.get('function')).replace ("['","").replace("']","")
						geneproduct[j]=str(gb_feature.qualifiers.get('product')).replace ("['","").replace("']","")
				gbfile.close()
				print "peptide_length", '\t', "downstream_gene", '\t', "downstream_gene_locus_tag", "\t", "downstream_gene_function",'\t',"downstream_gene_product","\t","blast_find_gene", '\t' "e_value","\t","query_length"
				blast_records = NCBIXML.parse(blast_file)
				for blast_record in blast_records:
					k=1
					for alignment in blast_record.alignments:
						for hsp in alignment.hsps:
							if hsp.expect < E_VALUE:
								downstream=blast_record.query
								acidlength=downstream.split()[0]
								gbkID=downstream.split()[4]
								locus_tag=downstream.split()[5]
								gbkID=int(gbkID)
								if 'hypothetical protein' not in alignment.title:
									if k==1:			
										title=alignment.title
										querylength=len(hsp.query)
										print acidlength,'\t',genename[gbkID],'\t',locus_tag,'\t',genefunction[gbkID],'\t',geneproduct[gbkID],'\t',title,'\t',hsp.expect,'\t',querylength
										k=k+1
				file.close()
				blast_file.close()

