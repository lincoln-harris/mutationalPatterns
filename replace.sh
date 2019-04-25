#!/bin/bash

for f in *; do
	#echo $f
	sed -i 's/#CHROM/#CHROM/g' $f
done