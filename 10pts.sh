#!/bin/bash  
# 10pts in the paper
dataPATH=/Users/guangzhao/Google-Drive/git/2018_exact_EHVI/data
methodPATH3=/Users/guangzhao/Google-Drive/git/2018_exact_EHVI/paper-grid-10d
methodPATH24=/Users/guangzhao/Google-Drive/git/2018_exact_EHVI/paper-wfg-CDD-AVL

cd $methodPATH3
make   #you can use make -f $methodPATH3/makefile, but that need absolute address in makefile
cd ..

cd $methodPATH24
make
cd ..

rm $dataPATH/time_grid.txt
rm $dataPATH/time_irs.txt
rm $dataPATH/time_wfg.txt

for value in 3d 4d 5d 6d 7d
do
    for k in {1..5}
    do
        echo "Begin running method grid" "$value"
        $methodPATH3/EHVI $dataPATH/ran.10pts.$value.10
    done
    for k in {1..10}
    do
        echo "Begin running method wfg & CDD" "$value"
        $methodPATH24/wfg1 $dataPATH/ran.10pts.$value.10
    done

done

for k in {1..10}
do
    echo "Begin running method wfg & CDD for 8d"
    $methodPATH24/wfg1 $dataPATH/ran.10pts.8d.10
done

echo "I am done running"
