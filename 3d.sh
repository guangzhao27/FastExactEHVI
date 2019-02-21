#!/bin/bash  
# 3d in the paper
dataPATH=/Users/guangzhao/Google-Drive/git/2018_exact_EHVI/data
timePATH =/Users/guangzhao/Google-Drive/git/2018_exact_EHVI/time
methodPATH12=/Users/guangzhao/Google-Drive/git/2018_exact_EHVI/paper-3d-grid-and-IRS
methodPATH34=/Users/guangzhao/Google-Drive/git/2018_exact_EHVI/paper-wfg-CDD-AVL
methodPATH5=/Users/guangzhao/Google-Drive/git/2018_exact_EHVI/paper-3d-KMAC

cd $methodPATH12
make   #you can use make -f $methodPATH3/makefile, but that need absolute address in makefile
cd ..

cd $methodPATH34
make
cd ..

cd $methodPATH5
make
cd ..

rm $timePATH/time_grid.txt
rm $timePATH/time_irs.txt
rm $timePATH/time_wfg.txt
rm $timePATH/time_cdd.txt
rm $timePATH/time_KMAC.txt
rm $timePATH/time_avl.txt

for value in 10pts 50pts 100pts 150pts 200pts 250pts 300pts
do
    for k in {1..5}
    do
        echo "Begin running method irs & grid" "$value"
        $methodPATH12/EHVI $dataPATH/ran.$value.3d.10

        echo "Begin running method wfg & cdd & irs" "$value"
       $methodPATH34/wfg1 $dataPATH/ran.$value.3d.10

        echo "Begin running method kmac" "$value"
        $methodPATH5/EHVI $dataPATH/ran.$value.3d.10
    done
done
echo "I am done running"
