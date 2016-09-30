#!/bin/sh

ssh -p 2222 mcastaldo@www.isislab.it

for((X=9;X<=24; X++))
do
ssh oxygen$X
done

exit 0
