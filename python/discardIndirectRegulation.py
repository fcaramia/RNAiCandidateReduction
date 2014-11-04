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
		
		if REG in db:
			db[REG].append(GENE)
		else:
			db[REG] = [GENE]

	return db

def readCandidates(candidate_file,candidates):
	
	candidates_obj = open(candidate_file)
	#Skip header
	next(candidates_obj)
	
	for line in candidates_obj:
		[GENE,SCORE] = line.rstrip('\n').split(',')
		
		candidates[GENE] = int(SCORE)

	return candidates

def discardCandidates(candidates,marks,reg_db,string_db, std_dev_filter):
	scores = []
	discarded = []
	for m in marks:
		for g in marks[m]:
			if candidates[m]>candidates[g]:
				scores.append(math.fabs(candidates[m]-candidates[g]))
	
	std = numpy.std(scores)

	for m in marks:
		for g in marks[m]:
			if candidates[m]>candidates[g]:
				if(math.fabs(candidates[m]-candidates[g])>=std_dev_filter*std):
					if m not in discarded: 
						discarded.append(m)

	return discarded	

def main():

	gene_db_dir = "../gene_db/"
	parser = optparse.OptionParser()
	parser.add_option('-t', '--transFac',action='store_true',dest="trans_fac",help="use TransFac DB")
	parser.add_option('-p', '--phosphoSite',action='store_true',dest="phospho_site",help="use Phosphosite DB")
	parser.add_option('-s', '--string',action='store_true',dest="string_db",help="use String DB")
	parser.add_option('-c', '--candidates',type='string',dest="candidate_file",help="candidates file")
	parser.add_option('-e', '--score', type='int',dest='score',help="minimum evidence score, default: 400" , default=400)
	parser.add_option('-d', '--direction', action='store_true',dest='direction',help="evidence of direction for interaction",default=False)
	parser.add_option('-x', '--standardDev',dest='std_dev',help="number of standard deviations to mark candidate for removal, default: 2",default=2.0, type='float')
	
	(options, args) = parser.parse_args()

	if not options.candidate_file:
		parser.error("candidate file is needed")
	
	#Read Candidates
	candidates = {}
	candidates = readCandidates(options.candidate_file,candidates)

	#Read DBs
	reg_db = {}
	marks = {}
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
		string_db = {}
		string_db = readStringActionDB(gene_db_dir+"actions_curated.tsv",string_db,options.score,options.direction)
		checkStringDb(candidates,string_db,marks)
	
	i = 0
	for r in marks:
		i += len(marks[r])
	print i, "candidates marked"

	if len(marks)>0:
		discarded = discardCandidates(candidates,marks,reg_db,string_db,options.std_dev)

	print "Discarded Genes:"
	for g in discarded:
		print g

	i = 0	
	print "Direct Regulators:"
	for c in candidates :
		if c not in discarded:
			i += 1
	print i 
	
if __name__ == '__main__':
	main()