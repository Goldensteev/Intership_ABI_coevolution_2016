#echo 'coucou'
#!/bin/bash -f 



if [ $# != 3 ]; then
	echo 'USAGE: '$0' <input alignement file (fst)> <input tree file (newick)> <output directory>'
	exit 1
fi



fiInAli=$1
fiInTree=$2
repOut=$3

ici=`pwd`


FASTML_EXE=/Users/steven/Desktop/Stage2016/FastML.v3.1/programs/fastml/fastml

if ! [ -e $FASTML_EXE ]; then
	echo "$0:: $FASTML_EXE not found"
	exit 1
fi

###############
if ! [ -e $fiInAli ]; then
	echo "$0:: $fiInAli not found"
	exit 1
fi
#Passer fiIn en path absolu
cd ${fiInAli%/*}
fiInAli=`pwd`"/"${fiInAli##*/}
cd $ici

################
if ! [ -e $fiInTree ]; then
	echo "$0:: $fiInTree not found"
	exit 1
fi
#Passer fiIn en path absolu
cd ${fiInTree%/*}
fiInTree=`pwd`"/"${fiInTree##*/}
cd $ici
##################

#Passer repOut en path absolu
if ! [ -d $repOut ]; then
	mkdir -p $repOut
fi
cd $repOut
repOut=`pwd`
cd $ici


cd $repOut

cmd="$FASTML_EXE -s $fiInAli -t $fiInTree -mw -qf -g -f "
echo "RUNNING "$cmd
$cmd &> $repOut/fastml.err
cd $ici





