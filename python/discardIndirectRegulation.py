import os
import optparse
import string
import sys
import threading
from subprocess import Popen, call
from os.path import  join

def readDB (db_file, db):
	
	db_obj = open(db_file, 'r+')
	
	#Skip header
	next(db_obj)
	

	for line in db_obj:
		[TF,GENE] = line.rstrip('\n').split(',')
		
		if TF in db:
			db[TF].append(GENE)
		else:
			db[TF] = [GENE]

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
	parser = optparse.OptionParser()
	parser.add_option('-d', '--dataBase',type='string',dest="data_base",help="database file")
	parser.add_option('-c', '--candidates',type='string',dest="candidate_file",help="candidates file")
	(options, args) = parser.parse_args()

	#Read DB
	db = {}
	db = readDB(options.data_base,db)
	

	#Read Candidates
	candidates = {}
	candidates = readCandidates(options.candidate_file,candidates)

	#print db
	#print candidates
	i = 0
	for c in candidates:
		if c in db:
			#print "Found TF " + c
			for g in db[c]:
				if g in c:
					i+=1
	print i		
	



if __name__ == '__main__':
	main()