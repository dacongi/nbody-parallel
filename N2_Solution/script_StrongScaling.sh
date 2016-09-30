
if [ ! -f generate_output.txt ]
then
        echo "file  does not exist.";
else
        rm generate_output.txt
fi

if [ ! -f stderror.csv ]
then
        echo "file  does not exist.";
else
        rm stderrorStrongScaling.csv
fi

/oxygen/compilers/mpi/openmpi/1.8.2.jdk.1.8/bin/mpicc $1 -o $1.out

for((X=2;X<=16; X*=2))
do
/oxygen/compilers/mpi/openmpi/1.8.2.jdk.1.8/bin/mpirun -np $X -host oxygen9,oxygen10,oxygen11,oxygen12,oxygen13,oxygen14,oxygen15,oxygen16,oxygen17,oxygen18,oxygen19,oxygen20,oxygen21,oxygen22,oxygen23,oxygen24 $1.out $2 10 0.1 1 g 2>> stderrorStrongScaling.csv 
done

exit 0
