#!/bin/sh

files=$(find $1 -name "*_036*pileUp.pu")

for i in $files
do
dirname=$(dirname $i)
head -n 3176 $i | tail -n 3096 | cut -f 1,2,3,4,5,6 > $i".trimmed.pu" &
done

