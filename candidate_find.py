#!/usr/bin/python
#Filename: candidate_find.py

from Bio import SeqIO
from Bio.Seq import Seq
import sys
import os

path=os.getcwd()
# motif of the overlap gene
seq1= Seq("TGATG")
seq2=Seq("ATGA")
seq3=Seq("TAATG")


def filter(begin, end, motif, strand,genestart,geneend,genestrand,genetype,locus_tag):
	if strand==1:
		for geneid,genestart in genestart.items():		
			if (genestart>(begin-1) and genestart< end):
				if genestrand[geneid]==1:
					if genetype[geneid]=="CDS":
						locus=str(locus_tag[geneid]).replace ("['","").replace("']","")
						print("%i:%i\t%s\t%i\t%i\t%s" \
						% (begin+1,end,motif,1,geneid,locus))
						break
	else:
		for geneid,geneend in geneend.items():					
			if (geneend>begin and geneend<end+1):
				if genestrand[geneid]==-1:
					if genetype[geneid]=="CDS":
						locus=str(locus_tag[geneid]).replace("['","").replace("']","")
						print("%i:%i\t%s\t%i\t%i\t%s" \
							% (begin+1,end,motif,-1,geneid,locus))
						break


def find_overlap_codon(seq,motif,genestart,geneend,genestrand,genetype,locus_tag):
		for i in range(0, len(seq)-1):
			if str(seq[i:i+len(motif)])==str(motif):
				filter(i,i+len(motif),motif,1,genestart,geneend,genestrand,genetype,locus_tag)
			if str(seq[i:i+len(motif)])==str(motif.reverse_complement()):

				filter(i,i+len(motif),motif,-1,genestart,geneend,genestrand,genetype,locus_tag)

for x in os.listdir(path):
	if ".gbk" in x:
		gb_file=x   #genbank file of the specie
		specie=x.split('.')[0]
		outputfilename=str(specie)+"_candidate.txt" 
		
		#the result file
		outfile=open(outputfilename,"w")
		infile=open(gb_file,"r")
		sys.stdout=outfile
		print "#Name", gb_file
		print "#location","\t","motif","\t","strand","\t","downstream_gbkid","\t","downstream_genelocus_tag"
		#hash				 
		genestart={}
		geneend={}
		genestrand={}
		locus_tag={}
		genetype={}
		for gb_record in SeqIO.parse(infile, "genbank") :
			number=len(gb_record.features)
			for j in range(0,number):
				gb_feature = gb_record.features[j]
				genestart[j] = gb_feature.location.nofuzzy_start
				geneend[j] = gb_feature.location.nofuzzy_end
				genestrand[j]=gb_feature.strand
				genetype[j]=gb_feature.type
				locus_tag[j]=gb_feature.qualifiers.get('locus_tag')
		find_overlap_codon(gb_record.seq,seq1,genestart,geneend,genestrand,genetype,locus_tag)
		find_overlap_codon(gb_record.seq,seq2,genestart,geneend,genestrand,genetype,locus_tag)
		find_overlap_codon(gb_record.seq,seq3,genestart,geneend,genestrand,genetype,locus_tag)
		outfile.close()
		infile.close()



