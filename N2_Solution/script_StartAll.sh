#!/bin/sh

for((X=9;X<=24; X++))
do
ssh oxygen$X
done

exit 0
