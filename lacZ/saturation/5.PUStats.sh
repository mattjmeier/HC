#!/bin/sh

files=$(find $1 -name "*trimmed.pu")

for i in $files
do
/utils/installs/cran/R_3_0_0/bin/Rscript  /home/mmeier/DNA/scriptsAndResults/lacZsaturation/bap/statsFromPU.R $i $i".stats" &
done
