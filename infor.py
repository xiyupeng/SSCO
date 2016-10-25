#!usr/bin/python
#filename: infor.py
from Bio import SeqIO
import sys
import os
#read the result xml file from blast

origin=sys.stdout
path=os.getcwd()

for x in os.listdir(path):
	if ("_annotated_upstream_gene" in x and "info" not in x):
		specie=x.split('_annotated')[0]
		infor_file=str(x+".info")
		file=open(infor_file,'w')
		sys.stdout=file
		readfile=open(x,'r')
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
				for line in readfile:
					if '#' in line:
						if 'NC' in line:
							print line
							print "downstream_gene", '\t', "downstream_gene_locus_tag", "\t", "downstream_gene_function",'\t',"downstream_gene_product","\t","upstream_gene", '\t', "upstream_gene_locus_tag", "\t", "upstream_gene_function",'\t',"upstream_gene_product"
					else:
						downgbkID=line.split()[3]
						downgbkID=int(downgbkID)
						downlocus_tag=line.split()[4]
						upgbkID=line.split()[5]
						upgbkID=int(upgbkID)
						uplocus_tag=line.split()[6]
						print genename[downgbkID], '\t', downlocus_tag, "\t", genefunction[downgbkID],'\t',geneproduct[downgbkID],"\t",genename[upgbkID], '\t', uplocus_tag, "\t", genefunction[upgbkID],'\t',geneproduct[upgbkID]
				readfile.close()
				file.close()
