import os
import optparse
import string
import sys
import threading
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
        self.gene_b_id = gene_d_id
        self.action = action
        self.mode = mode
        self.a_is_acting = a_is_acting
        self.score = score
        self.source = source

#read String db by certain filters    
def readStrinActionDB(db_file,db,min_score,direction,mode,action):
	db_obj = open(db_file, 'r+')
	#Skip header
	next(db_obj)
	for line in db_obj:
		[A,B,mode,action,a_is_acting,score,source,source2] = line.rstrip('\n').split(',')
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
		
		if GENE in candidates:
			candidates[GENE].append(SCORE)
		else:
			candidates[GENE] = [SCORE]

	return candidates

def main():

	gene_db_dir = "../gene_db"
	parser = optparse.OptionParser()
	parser.add_option('-t', '--transFac',action='store_true',dest="trans_fac",help="use TransFac DB")
	parser.add_option('-p', '--phosphoSite',action='store_true',dest="phospho_site",help="use Phosphosite DB")
	parser.add_option('-s', '--string',action='store_true',dest="string_db",help="use String DB")
	parser.add_option('-c', '--candidates',type='string',dest="candidate_file",help="candidates file")
	(options, args) = parser.parse_args()


	if not options.candidate_file:
		parser.error("candidate file is needed")
	#Read DB
	reg_db = {}
	
	if options.trans_fac:
		reg_db = readRegDB(options.trans_fac,reg_db)

	if options.phospho_site:
		reg_db = readRegDB(options.phospho_site,reg_db)

	#Read Candidates
	candidates = {}
	candidates = readCandidates(options.candidate_file,candidates)

	
	i = 0
	for c in candidates:
		if c in reg_db:
			#print "Found TF " + c
			for g in reg_db[c]:
				if g in c:
					i+=1
	print i		
	
if __name__ == '__main__':
	main()