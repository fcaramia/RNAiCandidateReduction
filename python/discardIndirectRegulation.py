import time
import os
import optparse
import string
import sys
import threading
import numpy
import math
from subprocess import Popen, call
from os.path import  join



class Interaction:

	gene_a_id = ""
	gene_b_id = ""
	action = ""
	source = ""
	mode = ""
	a_is_acting = False
	score = 0
	def __init__(self, gene_a_id,gene_b_id,mode,action,a_is_acting,score,source):
		self.gene_a_id = gene_a_id
		self.gene_b_id = gene_b_id
		self.action = action
		self.mode = mode
		self.a_is_acting = a_is_acting
		self.score = score
		self.source = source


def buildGraph(graph,reg_db,string_db):
	
	for gene in reg_db:
		if gene in graph:
			for i in reg_db[gene]:
				if i not in graph[gene]:
					graph[gene].append(i)
				else:
					continue
		else:
			graph[gene] = reg_db[gene]

	for gene in string_db:
		if gene in graph:
			for i in string_db[gene]:
				if i.gene_b_id not in graph[gene]:
					graph[gene].append(i.gene_b_id)
		else:
			graph[gene] = list(set([x.gene_b_id for x in string_db[gene]]))

	return graph


def havePath(graph,start,end, max_depth,path=[],visited=[],depth=0):
	
	
	if depth > max_depth:
		return None
	path = path + [start]
	if start == end :
		return path
	if start not in graph:
		return None
	visited.append(start)
	for node in graph[start]:
		if node in visited:
			continue
		ret = havePath(graph,node,end,max_depth,path,visited,depth+1)
		if ret: return ret

	return None
			

def checkGraph(candidates,graph,marks,max_depth):

	part =  (len(candidates)**2)/100
	count = 0
	for i in range(len(candidates)):
		for j in range(i+1, len(candidates)):
			count = count + 1
			
			if count%part == 0:
				print count/part,"%"

			c1 = candidates[i]
			c2 = candidates[j]
			#print "Evaluate candidates: ",c1,c2
			if c1 in marks and c2 in marks[c1]:
				continue 
			if c2 in marks and c1 in marks[c2]:
				continue	
			

			path = havePath(graph,c1,c2,max_depth,visited=[])
			if path :
				
				if c1 in marks:
					marks[c1].append(c2)
				else:
					marks[c1] = [c2]

			else:
				path = havePath(graph,c2,c1,max_depth,visited=[])
			
				if path:
					
					if c2 in marks:
						marks[c2].append(c1)
					else:
						marks[c2] = [c1]

			#if path:
			#	print path
	return marks
#Check if regulation interaction is present in candidates
#Check for repeats and self-regulation
def checkRegDb(candidates,db,marks):
	for c in candidates:
		if c in db:
			for g in db[c]:
				if g!=c and g in candidates:
					if c in marks and g not in marks[c]:
						marks[c].append(g)
					else:
						marks[c]=[g]
	return marks

def checkStringDb(candidates,db,marks):	
	for c in candidates:
		if c in db:
			for inter in db[c]:
				if inter.gene_b_id!=c and inter.gene_b_id in candidates:
					if c in marks and inter.gene_b_id not in marks[c]:
						marks[c].append(inter.gene_b_id)
					else:
						marks[c]=[inter.gene_b_id]
	return marks
#read String db by certain filters    
def readStringActionDB(db_file,db,min_score,direction,mode=None,action=None):
	db_obj = open(db_file, 'r+')
	#Skip header
	next(db_obj)
	for line in db_obj:
		[A,B,mode,action,a_is_acting,score,source,source2] = line.rstrip('\n').split('\t')
		#Apply filters
		if min_score > int(score):
			continue
		if A == ' ' or B == ' ':
			continue
		if A == '' or B == '':
			continue
		if A == "NA" or B == "NA":
			continue
		if direction and a_is_acting == '0':
			continue
		inter = Interaction(A,B,mode,action,a_is_acting,score,source)
		
		if A in db :
			db[A].append(inter)
		else:
			db[A] = [inter]

	return db


#Same function for transfac and phosphosite
def readRegDB (db_file, db):
	
	db_obj = open(db_file, 'r+')
	#Skip header
	next(db_obj)
	for line in db_obj:
		[REG,GENE] = line.rstrip('\n').split(',')
		if REG == ' ' or GENE == ' ':
			continue
		if GENE == '' or REG == '':
			continue
		if GENE == "NA" or REG == "NA":
			continue
		if REG in db:
			if GENE not in db[REG]:
				db[REG].append(GENE)
		else:
			db[REG] = [GENE]

	return db

def readCandidates(candidate_file,candidates):
	
	candidates_obj = open(candidate_file)
	#Skip header
	next(candidates_obj)
	scores = []
	for line in candidates_obj:

		[GENE,SCORE] = line.rstrip('\n').split(',')
		
		candidates[GENE] = float(SCORE)
		
	return candidates

def discardCandidates(candidates,marks,reg_db,string_db, std_dev_filter, zscore):
	scores = []
	discarded = []
	for m in marks:
		for g in marks[m]:
			if candidates[m]>candidates[g]:
				scores.append(math.fabs(candidates[m]-candidates[g]))
			

	std = numpy.std(scores)

	for m in marks:
		for g in marks[m]:
			if candidates[m]<candidates[g]:
				if(math.fabs(candidates[m]-candidates[g])>=std_dev_filter*std):
					if m not in discarded and candidates[m]<=zscore: 
						discarded.append(m)
			

	return discarded	

def main():

	gene_db_dir = "../gene_db/"
	parser = optparse.OptionParser()
	parser.add_option('-t', '--transFac',action='store_true',dest="trans_fac",help="use TransFac DB")
	parser.add_option('-p', '--phosphoSite',action='store_true',dest="phospho_site",help="use Phosphosite DB")
	parser.add_option('-s', '--string',action='store_true',dest="string_db",help="use String DB")
	parser.add_option('-c', '--candidates',type='string',dest="candidate_file",help="candidates file")
	parser.add_option('-r', '--score', type='int',dest='score',help="minimum evidence score, default: 400" , default=400)
	parser.add_option('-d', '--direction', action='store_true',dest='direction',help="evidence of direction for interaction",default=False)
	parser.add_option('-x', '--standardDev',dest='std_dev',help="number of standard deviations to mark candidate for removal, default: 2",default=2.0, type='float')
	parser.add_option('-z', '--maxZScore',dest='zscore',help="max zscore of discarded gene",default=3.0, type='float')
	parser.add_option('-o', '--outfile',dest='out',help="output file",default='out.txt', type='string')
	parser.add_option('-e', '--expand',dest='expand',help="create graph to expand regulation assumption",action="store_true",default=False)
	parser.add_option('-m', '--maxExpand',dest='max_expand',help="max size of path, default: 5",default=5, type=int)


	(options, args) = parser.parse_args()

	if not options.candidate_file:
		parser.error("candidate file is needed")
	
	
	#Read Candidates
	candidates = {}
	candidates = readCandidates(options.candidate_file,candidates)
	
	scores = []
	for c in candidates:
		scores.append(candidates[c])

	sign = numpy.mean(scores)

	#Validate candidates	 
	for c in candidates:
		if sign > 0 and candidates[c]<0:
			parser.error("z-scores should have same sign")

		if sign < 0 and candidates[c]>0:
			parser.error("z-scores should have same sign")	


	##Change direction of regulation if negative
	if sign <0:
		for c in candidates:
			candidates[c]*=-1


	#Read DBs
	graph = {}
	reg_db = {}
	marks = {}
	string_db = {}
	if options.trans_fac:
		print "Loading TransFac"
		reg_db = readRegDB(gene_db_dir+"transfac_interactions.csv",reg_db)

	if options.phospho_site:
		print "Loading Phosphosite"
		reg_db = readRegDB(gene_db_dir+"kinase_curated_db.csv",reg_db)

	if options.phospho_site or options.trans_fac:
		marks = checkRegDb(candidates,reg_db,marks)

	if options.string_db:
		print "Loading String Actions"
		string_db = readStringActionDB(gene_db_dir+"actions_curated.tsv",string_db,options.score,options.direction)
		marks = checkStringDb(candidates,string_db,marks)
	
	if options.expand:
		if options.direction is False:
			parser.error("To expand, direction must be used")
		print "Expanding Search"
		print "Building Graph"
		graph = buildGraph(graph, reg_db, string_db)
		print "Checking for paths"
		marks = checkGraph(candidates.keys(), graph, marks,options.max_expand)


	i = 0
	for r in marks:
		i += len(marks[r])
	print i, "interactions marked"

	f = open(options.out, 'w+')
	f.write("Input file: "+options.candidate_file+'\n')	
	f.write("DB Interactions marked: "+str(i)+'\n')
	if len(marks)>0:
		discarded = discardCandidates(candidates,marks,reg_db,string_db,options.std_dev,options.zscore)
		f.write("Discarded Genes:\n")
		for r in discarded:
			f.write(r+' '+str(candidates[r])+ "\n")
		
		f.write("\n\nDiscarded Genes detailed:\n")
		for r in discarded:
			f.write(r+' '+str(candidates[r])+ " regulates: \n")
			for g in marks[r]:
				f.write("\t"+g+" "+str(candidates[g])+'\n')
		print "Discarded Genes: ",len(discarded)
		
	i = 0	
	print "Direct Regulators:"
	f.write("\n\nGenes not dicarded:\n")
	for c in candidates :
		if c not in discarded:
			f.write(c+' '+str(candidates[c])+'\n')
			i += 1
	print i 
	
if __name__ == '__main__':
	main()