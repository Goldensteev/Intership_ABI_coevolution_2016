for i in `cat listeAli`; do sed "s/MONALI/$i/g" ./LoadLeveler_FASTML.txt > 1.submit ; llsubmit 1.submit; done

for i in `cat listTrees`; do sed "s/MONTREE/$i/g" ./LoadLeveler_epics.txt > 2.submit ; llsubmit 2.submit; done

for i in `cat listallCATH`; do sed "s/MONDOMPDB/$i/g" ./distcalc.txt > 3.submit;

