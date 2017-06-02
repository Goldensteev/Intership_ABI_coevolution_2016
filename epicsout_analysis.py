# -*- coding: utf-8 
"""
Created on Mon Mar 21 10:22:16 2016

@author: steven
"""

import os
import re
import numpy as np
#np.set_printoptions(threshold=np.nan)


#function parses the outfile off epics-dist
#INPUT - the epics outfile 
#OUTPUT - a list done of the lines of the outfile      
def parse_epicsout(epicsout):
    f=open(epicsout,'r')
    out=f.readlines()
    f.close()
    return out



  
#function make matrix of the position and nb of chronologies and/or coocurences and the pvalues
#INPUT - the list output by parse_epicsout()
#OUTPUT - the matrix with nb of chronologies and/or coocurences and one with the pvalues 
def get_events(epicsout):
    seqlen=len(re.findall('(?<=\[)(.*?)(?=\])',epicsout[0])[0].split('/'))
    mat_nbevents=np.ones((seqlen,seqlen))
    mat_pval=np.ones((seqlen,seqlen))
    for h in epicsout[1:]:
        if 'Root' not in h and '|' not in h:
            k=re.findall(r'\d+',h)
            pos1=int(k[0])
            pos2=int(k[1])
            ch=int(k[2])
            if len(k)==5 :
                pval=float(k[3]+'.'+k[4])
            else :
                pval=1
            mat_nbevents[pos1][pos2]=ch
            mat_pval[pos1][pos2]=pval
    return mat_nbevents,mat_pval
    

####MAIN####
    
def do_all_files(RESpath,OUTpath):
    listfiles=os.listdir(RESpath)
    for i in listfiles:
        t=i.split('_')
        print t
        j=parse_epicsout(RESpath+i)
        l=get_events(j)
        np.savetxt(OUTpath+t[0]+'_mat_'+t[1],l[0],fmt='%.f')
        np.savetxt(OUTpath+t[0]+'_matpval_'+t[1],l[1],fmt='%.6f')


do_all_files('/Users/steven/Desktop/Stage2016/SpecificFams/new/epicsRes/raxml/','/Users/steven/Desktop/Stage2016/SpecificFams/new/epicsMat/raxml/')