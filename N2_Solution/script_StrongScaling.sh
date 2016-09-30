
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

mpicc $1 -o $1.out

for((X=2;X<=16; X*=2))
do
mpirun -np $X $1.out $2 8 1 1 g 2>> stderror.csv
done

exit 0
