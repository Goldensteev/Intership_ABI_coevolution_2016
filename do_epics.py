# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 13:54:31 2016

@author: steven
"""

import os


listfiles=os.listdir('/Users/steven/Desktop/Stage2016/SpecificFams/new/fastmlruns/raxml/')[1:]
print listfiles
for i in listfiles :
    j=i.split('_')[0]
    #os.system('python events.py /Users/steven/Desktop/Stage2016/SpecificFams/new/fastmlruns/raxml/'+i+'/ /Users/steven/Desktop/Stage2016/SpecificFams/new/epicsTrees/raxml/')
    os.system('/Users/steven/Desktop/Stage2016/Diversity/libtree/epics -I -a < /Users/steven/Desktop/Stage2016/SpecificFams/new/epicsTrees/raxml/'+j+'_epicsTree'+' >>/Users/steven/Desktop/Stage2016/SpecificFams/new/epicsRES/raxml/'+j+'_cooc')