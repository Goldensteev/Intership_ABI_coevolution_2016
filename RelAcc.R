maxsolvent=read.table('/Users/steven/Desktop/Stage2016/SpecificFams/new/analyses/maxsolvent.txt')
datatest=cleandata2
for(i in 1:length(datatest$SEQ)){
  datatest$RelAccAA1[i]=100*(datatest$enfouiAA1[i]/(maxsolvent[maxsolvent$V1==levels(droplevels(datatest$AA1[i])),]$V2))
  datatest$RelAccAA2[i]=100*(datatest$enfouiAA1.1[i]/(maxsolvent[maxsolvent$V1==levels(droplevels(datatest$AA2[i])),]$V2))
}
write.csv(datatest,file = '/Users/steven/Desktop/Stage2016/SpecificFams/new/analyses/relacc001/3.40.980.10_relacc_clean.csv')



dall1=read.table(file = '/Users/steven/Desktop/Stage2016/SpecificFams/new/ExempleFamilles/OK/2.40.180.10/AllRes.txt',sep=' ')
dall1b=dall1[1:4]
m1all=matrix(0,415,415)
for(i in 1:length(dall1b$V1)){
      
}