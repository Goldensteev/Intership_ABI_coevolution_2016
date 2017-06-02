library(lattice)
library(car)

minpositive = function(x) min(x[x > 0])

data1=read.table('/Users/steven/Desktop/Stage2016/SpecificFams/new/RES/raxml/ 2.40.180.10',header = TRUE)
data2=read.table('/Users/steven/Desktop/Stage2016/SpecificFams/new/RES/raxml/ 3.10.200.10',header = TRUE)
#data3=read.table('/Users/steven/Desktop/Stage2016/SpecificFams/new/RES/raxml/ 3.40.50.410',header = TRUE)
data4=read.table('/Users/steven/Desktop/Stage2016/SpecificFams/new/RES/raxml/ 3.40.800.10',header = TRUE)
data5=read.table('/Users/steven/Desktop/Stage2016/SpecificFams/new/RES/raxml/ 3.40.980.10',header = TRUE)


mat1=read.table('/Users/steven/Desktop/Stage2016/SpecificFams/new/epicsMat/raxml/2.40.180.10_matpval_cooc')
mat2=read.table('/Users/steven/Desktop/Stage2016/SpecificFams/new/epicsMat/raxml/3.10.200.10_matpval_cooc')
#mat3=read.table('/Users/steven/Desktop/Stage2016/SpecificFams/new/epicsMat/cleaned/3.40.50.410_matpval_cooc')
mat4=read.table('/Users/steven/Desktop/Stage2016/SpecificFams/new/epicsMat/raxml/3.40.800.10_matpval_cooc')
mat5=read.table('/Users/steven/Desktop/Stage2016/SpecificFams/new/epicsMat/raxml/3.40.980.10_matpval_cooc')

  
m1=as.matrix(mat1)
m2=as.matrix(mat2)
#m3=as.matrix(mat3)
m4=as.matrix(mat4)
m5=as.matrix(mat5)

t1=table(m1)
t2=table(m2)
#t3=table(m3)
t4=table(m4)
t5=table(m5)


plot(t1,main='2.40.180.10')
plot(t2,main='3.10.200.10')
plot(t3,main='3.40.50.410')
plot(t4,main='3.40.800.10')
plot(t5,main='3.40.980.10')
minpositive(data1[11])

t1struct=table(data1[5])+table(data1[9])
t2struct=table(data2[5])+table(data2[9])
t3struct=table(data3[5])+table(data3[9])
t4teststruct=table(datatest[5])+table(datatest[9])
sum(t1struct)


cleandata1=data1
cleandata1<-cleandata1[!(cleandata1$distpair==-1),]
cleandata1<-cleandata1[!(cleandata1$AA1=='-'),]
write.csv(cleandata1,file = '/Users/steven/Desktop/Stage2016/SpecificFams/new/analyses/raxml/2.40.180.10_clean.csv')

cleandata2=data2
cleandata2<-cleandata2[!(cleandata2$distpair==-1),]
cleandata2<-cleandata2[!(cleandata2$AA1=='-'),]
write.csv(cleandata2,file = '/Users/steven/Desktop/Stage2016/SpecificFams/new/analyses/raxml/3.10.200.10_clean.csv')

cleandata3=data3
cleandata3<-cleandata3[!(cleandata3$distpair==-1),]
cleandata3<-cleandata3[!(cleandata3$AA1=='-'),]
write.csv(cleandata3,file = '/Users/steven/Desktop/Stage2016/SpecificFams/new/analyses/raxml/3.40.50.410_clean.csv')

cleandata4=data4
cleandata4<-cleandata4[!(cleandata4$distpair==-1),]
cleandata4<-cleandata4[!(cleandata4$AA1=='-'),]
write.csv(cleandata4,file = '/Users/steven/Desktop/Stage2016/SpecificFams/new/analyses/raxml/3.40.800.10_clean.csv')

cleandata5=data5
cleandata5<-cleandata5[!(cleandata5$distpair==-1),]
cleandata5<-cleandata5[!(cleandata5$AA1=='-'),]
write.csv(cleandata5,file = '/Users/steven/Desktop/Stage2016/SpecificFams/new/analyses/raxml/3.40.980.10_clean.csv')

#rank
v1=as.vector(m1)
r1=rank(v1[v1<1])
plot(r1,v1[v1<1],main='Plot of Pvalue against Rank 2.40.180.10',ylab = 'Pvalue',xlab='rank')

v2=as.vector(m2)
r2=rank(v2[v2<1])
plot(r2,v2[v2<1],main='Plot of Pvalue against Rank 3.10.200.10',ylab = 'Pvalue',xlab='rank')

v4=as.vector(m4)
r4=rank(v4[v4<1])
plot(r4,v4[v4<1],main='Plot of Pvalue against Rank 3.40.800.10',ylab = 'Pvalue',xlab='rank')

v5=as.vector(m5)
r5=rank(v5[v5<1])
plot(r5,v5[v5<1],main='Plot of Pvalue against Rank 3.40.980.10',ylab = 'Pvalue',xlab='rank')


#tree reading & stuff
treetest=read.tree(file = '/Users/steven/Desktop/Stage2016/SpecificFams/new/newickTrees/2.40.180.10_tree.newick.txt')
rootedtree=root(treetest,outgroup = 41)
write.tree(rootedtree, file = '/Users/steven/Desktop/Stage2016/SpecificFams/new/treetest')
