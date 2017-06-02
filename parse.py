# -*- coding: utf-8 -*-
"""
Created on Fri May 20 13:30:24 2016

@author: steven
"""

import os 
from Bio import SeqIO


seq='----MFSEPIGVSVDGYPGIAIVVKNPPEFK-SIELRDEAIKFDAKGA----MILTVKDN--------EIVFPEDF----------RPLKEMYPEVAKKI----------------------------------------VDYEDGDAVIITWAETPAKALKSAIHVAYILKKEEI---------'
print seq

def parse_endres(resfile,Acctable):
    f=open(resfile,'r')
    res=f.read().splitlines()
    f.close()
    for i in res[1:] :
        j=i.split()
        if j==[] :
            continue
        ###Get_min###
        #getminmax_row(j,10)
        RelAcc1,RelAcc2=do_Relative_accessible_surface_area(j)


def getminmax_row(j,col):
    mindist=None
    minrow=None
    maxdist=None
    maxrow=None
    if float(j[col])!=-1 and float(j[col])<mindist or mindist==None :
        mindist=float(j[col])
        minrow=j
    if float(j[col])!=-1 and float(j[col])>maxdist or maxdist==None :
        maxdist=float(j[col])
        maxrow=j    
    return mindist,minrow,maxdist,maxrow



def get_acctable():
    f=open('/Users/steven/Desktop/Stage2016/maximum possible solvent accessible surface area.txt')
    Acc=f.read().splitlines()
    f.close()
    Acctable={}
    for i in Acc :
        #print i.split()
        Acctable[i.split()[0]]=i.split()[1]
    return Acctable    

    
def do_Relative_accessible_surface_area(j,Acctable):
    RelAcc1=0
    RelAcc2=0
    RelAcc1=100*j[5]/Acctable[j[2]]
    RelAcc2=100*j[9]/Acctable[j[6]]
    return RelAcc1,RelAcc2
    
    
Acctable=get_acctable()    
parse_endres('/Users/steven/Desktop/Stage2016/SpecificFams/Gyrase_S100.3.30.1360.30/ S100.3.30.1360.30',Acctable)        