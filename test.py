# -*- coding: utf-8 -*-
"""
Created on Mon May  9 13:17:21 2016

@author: steven
"""


import numpy as np

a=[[5,5,5,5,5],[5,5,5,5],[5,5,5],[5,5],[5]]
b='-AB-BC-AA'
listgaps=[]

def make_fullmat(a):
    for i in range(len(a)):
        for j in range(i+1):
            a[i].insert(0,0)
    a.insert(len(a),[0]*(len(a)+1))
    return a

def make_symetricalmat(a):
    for i in range(len(a)):
        for j in range(len(a)):
            a[j][i]=a[i][j]
    return a
    

def reformat_mat(oldmat,seq):
    nbgaps=0
    oldmat=np.asarray(oldmat)
    newmat=np.zeros((len(seq),len(seq)))
    for i in range(len(seq)):
        if seq[i]=='-':
            newmat[i]=-1
            newmat[:,i]=-1
            nbgaps+=1
        else :
            nbgapsbis=0
            for j in range(len(seq)) :
                if seq[j]=='-':
                    newmat[i][j]=-1
                    nbgapsbis+=1
                else :
                    newmat[i][j]=oldmat[i-nbgaps][j-nbgapsbis]
        newmat[i][i]=0
    return newmat
 


        
def getnbGAPS(seq,posAA):
    nbgaps=0
    for pos in range(len(seq)) :
        if seq[pos]=='-' :
            nbgaps+=1
        if pos >= posAA :
            return nbgaps



def do_len_trimming():
    f=open('/Users/steven/Desktop/Stage2016/listsmall','r')
    listall=f.read().splitlines()
    f.close()
    listsmall=[]
    print len(listall)
    for i in range(len(listall)):
        listall[i]=listall[i].split()
        print listall[i][0]
        if int(listall[i][0]) <= 100 :
            listsmall.append(listall[i])
    return listsmall


print do_len_trimming()
#print getnbGAPS('AAA--A',5)
#    
#a=make_fullmat(a)    
#a=make_symetricalmat(a)
#print a
#a=reformat_mat(a,b)
#print a
#
##print os.listdir('/Users/steven/Desktop/Stage2016/DATA/CATH_DATA/SequenceBySuperfamily/COMBS_CONSENSUS_TRIM/NR/MSA_MA/')[-1]
#    
##==============================================================================
## for i in range(len(b)):
##     if b[i]=='-':
#         a.insert(i,[-1]*len(b))
#         listgaps.append(i)
# 
# for i in range(len(b)):
#     if b[i]=='-':
#         for j in range(len(b)):
#             if j not in listgaps :
#                 a[j].insert(i,-1)
#         nbgaps+=1
# print a
#==============================================================================

#print np.asmatrix(a)
#
#print symmetrize(np.asmatrix(a))
#
#def symmetrize(a):
#    return a + a.T - numpy.diag(a.diagonal())
