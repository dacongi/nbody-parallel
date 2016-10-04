#!/bin/sh

if [ -d StrongScaling_CSV ]
        then
        rm -rf StrongScaling_CSV
        mkdir StrongScaling_CSV
        else
        mkdir StrongScaling_CSV
fi

/oxygen/compilers/mpi/openmpi/1.8.2.jdk.1.8/bin/mpicc $1 -o $1.out tree.c -lm

for((Y=50000; Y<=150000; Y+=50000))
do

        echo "particles,rank,time" >> stderrorStrongScaling_$Y.csv

        for((X=2;X<=16; X=X*2))
        do
        
                if [ -f generate_output.txt ]
                        then
                        rm generate_output.txt
                fi
        
        /oxygen/compilers/mpi/openmpi/1.8.2.jdk.1.8/bin/mpirun -np $X -host oxygen9,oxygen10,oxygen11,oxygen12,oxygen13,oxygen14,oxygen15,oxygen16,oxygen17,oxygen19,oxygen20,oxygen21,oxygen22,oxygen23,oxygen24 $1.out $Y 10 0.1 1 g 2>> stderrorStrongScaling_$Y.csv 
        
        done

        mv stderrorStrongScaling_$Y.csv StrongScaling_CSV/
done

exit 0
