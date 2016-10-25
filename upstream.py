#!/usr/bin/python
#Filename: upstream.py

#use the candidate.txt (result of the candidate_find.py) to find out the possible ORF
#we cut the result into two part: one is for the upstream gene that is annotated; the other is for the upstream that is unknown
#we only consider cds gene
from Bio import SeqIO
from Bio.Seq import Seq
import sys
import os

path=os.getcwd()
#if the gene is not annotated
def unknowngene(gb_file,begin,end,strand,motif,gbkID,locus,genestart,geneend,genestrand,genetype,record):
	k=1
	while k>0:
		if strand==1:
			newID=gbkID-k
		else:
			newID=gbkID+k
		if (newID==-1 or newID>=len(genetype)):
			return 0
		if (genetype[newID]=="CDS" and genestrand[newID]==strand):
			lastid=newID
			break
		k=k+1

	if strand==1:
		bound=int((begin-geneend[lastid]-1)/3)
		if 'ATGA' in motif:
			biggesti=0
			for i in range(1,bound):                                  #find the longest ORF
				codoncheck=str(record[begin-3*i:begin+3-3*i])
				if (codoncheck=='ATG'or codoncheck=='CTG'):
					biggesti=i
				if (codoncheck=='TAA'or codoncheck=='TAG' or codoncheck=='TGA'):
					break
			if biggesti>29:
				print '>',biggesti,location,motif,strand,gbkID,locus
				print str(record[begin-3*biggesti:begin+3])
		else:
			biggesti=0
			for i in range(1,bound):
				codoncheck=str(record[begin-1-3*i:begin+2-3*i])
				if (codoncheck=='ATG'or codoncheck=='CTG'):
					biggesti=i
				if (codoncheck=='TAA'or codoncheck=='TAG' or codoncheck=='TGA'): 
					break
			if biggesti>29:
				print '>',biggesti,location,motif,strand,gbkID,locus
				print str(record[begin-3*biggesti-1:begin+2])

	if strand==-1:
		bound=int((genestart[lastid]-end)/3)
		if 'ATGA' in motif:
			biggesti=0
			for i in range(1,bound):                                  #find the longest ORF
				codoncheck=str(record[end-3+3*i-1:end+3*i-1])
				if (codoncheck=='CAT'or codoncheck=='CAG'):
					biggesti=i
				if (codoncheck=='TTA'or codoncheck=='CTA' or codoncheck=='TCA'):
					break
			if biggesti>29:
				print '>',biggesti,location,motif,strand,gbkID,locus
				print str(record[end-4:end+3*biggesti-1])
		else:
			biggesti=0
			for i in range(1,bound):
				codoncheck=str(record[end-3+3*i:end+3*i])
				if (codoncheck=='CAT'or codoncheck=='CAG'):
					biggesti=i
				if (codoncheck=='TTA'or codoncheck=='CTA' or codoncheck=='TCA'): 
					break
			if biggesti>29:
				print '>',biggesti,location,motif,strand,gbkID,locus
				print str(record[end-3:end+3*biggesti])

#if the gene is annotated

def gene(gb_file,begin,end,strand,motif,gbkID,locus,genestart,geneend,genestrand,locus_tag,genetype):
	
	k=1
	while k>0:
		if strand==1:
			newID=gbkID-k
		else:
			newID=gbkID+k
		if (newID==-1 or newID>=len(genetype)):
			return 0
		if (genetype[newID]=="CDS" and genestrand[newID]==strand):
			lastid=newID
			break
		k=k+1
	if strand==1:
		if (geneend[lastid]>begin and geneend[lastid]<end):
			uplocus=str(locus_tag[lastid]).replace("['","").replace("']","")
			location=str(begin)+":"+str(end)
			print location,'\t',motif,'\t',strand,'\t',gbkID,'\t',locus,'\t',lastid,'\t',uplocus
			return 1
		else:
			return 0
	if strand==-1:
		if (genestart[lastid]>begin-1 and genestart[lastid]<end-1):
			uplocus=str(locus_tag[lastid]).replace("['","").replace("']","")
			location=str(begin)+":"+str(end)
			print location,'\t',motif,'\t',strand,'\t',gbkID,'\t',locus,'\t',lastid,'\t',uplocus
			return 1
		else:
			return 0

origin=sys.stdout
for x in os.listdir(path):
	if "_candidate.txt" in x:
		specie=x.split('_candidate')[0]
		anno_file=specie+'_annotated_upstream_gene'
		anno=open(anno_file,'w')
		unkn_file=specie+'_unknown_upstream_gene.fa'
		unknown=open(unkn_file,'w')
		candidate=open(x,'r')
		gb_file=str(specie+".gbk")
		for y in os.listdir(path):
			if gb_file in y:
				genestart={}
				geneend={}
				genestrand={}
				locus_tag={}
				genetype={}
				gbfile=open(gb_file,"r")
				for gb_record in SeqIO.parse(gbfile, "genbank") :
					record=gb_record.seq
					number=len(gb_record.features)
					for j in range(0,number):
						gb_feature = gb_record.features[j]
						genestart[j] = gb_feature.location.nofuzzy_start
						geneend[j] = gb_feature.location.nofuzzy_end
						genestrand[j]=gb_feature.strand
						genetype[j]=gb_feature.type
						locus_tag[j]=gb_feature.qualifiers.get('locus_tag')
				gbfile.close()
				for line in candidate:
					line=line.strip('\n')
					if '#' in line:
						if 'NC' in line:
							sys.stdout=anno
							print line,"annotated_upstream_gene"
							seqprint=("#location","motif","strand","downstream_gbkID","downstreamLocus_tag","Upstream_gbkID","upstreamLocus_tag")
							print "\t".join(seqprint)
					else:
						location,motif,strand,gbkID,locus=line.split('\t')
						begin,end=location.split(':')
						begin=int(begin)
						end=int(end)
						gbkID=int(gbkID)
						strand=int(strand)
						sys.stdout=anno
						find=gene(gb_file,begin,end,strand,motif,gbkID,locus,genestart,geneend,genestrand,locus_tag,genetype)
						if find == 0:
							sys.stdout=unknown
							re=unknowngene(gb_file,begin,end,strand,motif,gbkID,locus,genestart,geneend,genestrand,genetype,record)
				candidate.close()
				anno.close()
				unknown.close()
