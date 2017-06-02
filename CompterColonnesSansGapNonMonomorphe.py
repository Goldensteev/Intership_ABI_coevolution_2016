## -*- coding: utf-8 -*-

#!/bin/python
import sys
import re
#import RendreAliIdentiqueToPDB
#import SeqAliIdentiques

def OpenAndReadFile(nomFi):
	try:
		f = open(nomFi, "r")
	except IOError:
		print "OpenAndReadFile:: Fichier <%s> introuvable, arret du programme"%(nomFi)
		sys.exit(1) 
	t=f.read()
	ll=t.splitlines()
	f.close()
	return ll


def readFastaMul(nomFi):
	#Lecture du contenu du fichier
	lines=OpenAndReadFile(nomFi)
	#Parsage du du contenu du fichier
	seq=""
	nom=""
	lesSeq={}
	for l in lines:
		if l[0]==">":
			if seq != "":
				#tmp=(nom,seq)
				lesSeq[prot]=seq
				#print "adding ",prot
			if l[0:5]==">cath":
				t2=l.split("_")
				prot=t2[4]
			else:
				prot=l[1:7]
				
			seq=""
		else:
			#Ne pas oublier d'enlever le retour chariot a la fin des lignes
			seq=seq+l[:-1]
	if seq != "":
		#tmp=(nom,seq)
		lesSeq[prot]=seq
	return lesSeq

# one_letter["SER"] will now return "S"
one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
'GLY':'G', 'PRO':'P', 'CYS':'C'}

# three_letter["S"] will now return "SER"
three_letter = dict([[v,k] for k,v in one_letter.items()])

three_letter ={'V':'VAL', 'I':'ILE', 'L':'LEU', 'E':'GLU', 'Q':'GLN', \
'D':'ASP', 'N':'ASN', 'H':'HIS', 'W':'TRP', 'F':'PHE', 'Y':'TYR',    \
'R':'ARG', 'K':'LYS', 'S':'SER', 'T':'THR', 'M':'MET', 'A':'ALA',    \
'G':'GLY', 'P':'PRO', 'C':'CYS'}

accessibility={"A":107.95, "C":134.28, "D":140.39, "E":172.25, "F":199.48, "G":80.10, "H":182.88, "I":175.12, "K":200.81, "L":178.63, "M":194.15, "N":143.94, "P":136.13, "Q":178.50, "R":238.76, "S":116.50, "T":139.27, "V":151.44, "W":249.36, "Y":212.76}
accessibility3={}
for i in accessibility:
	accessibility3[three_letter[i]]=accessibility[i]


##########################################



#>cath|4_0_0|1cf9A01/51-436
#
######Lecture fichier ali mul:
lesSeqAli=readFastaMul(sys.argv[1])
lesProt=lesSeqAli.keys()
lgAli=len(lesSeqAli[lesProt[0]])
nbProt=len(lesProt)

colonneOK=[True]*len(lesSeqAli[lesProt[0]])
#Chaque protein a une liste de booleen ou pour chauqe aa le booleen dit s'il est dans un colonne sans gap ou non
posOK={}
nbColOK=0
for p in lesProt:
	posOK[p]=[]
for i in range(lgAli):
	#D'abord regarder s'il y a un gap dans la colonne
	#Et regarder si la colonne est monomorphe
	j=0
	aaDifferent=False
	aa=lesSeqAli[lesProt[0]][i]
	while colonneOK[i] and j<nbProt:
		if lesSeqAli[lesProt[j]][i]=='-':
			colonneOK[i]=False
		if lesSeqAli[lesProt[j]][i]!=aa: #On met à True des qu'on trouve un aa different
			aaDifferent=True
		j+=1
	if colonneOK[i]: # si pas de gap, met False si colonne monomorqphe
		colonneOK[i]=aaDifferent
	if colonneOK[i]: # si pas de gap, et  si pas  monomorqphe
		nbColOK+=1
		
print sys.argv[1],"nb_colonnes_OK", nbColOK, "nb_Paires", nbColOK*(nbColOK-1)/2