
# Make BED file with positions of all violations
cat $violations | grep -v "^#" | awk ' { print $1"\t"$2-1"\t"$2 } ' > /tmp/$child.bed


echo "Running GATK DepthOfCoverage..."
while true
do
	java -jar $GATK3446    -T DepthOfCoverage -R $mm10min    -o /tmp/$child.sire.allele_depth    -I /tmp/$sire.list --printBaseCounts    -L /tmp/$child.bed
	pid1=$!
	java -jar $GATK3446    -T DepthOfCoverage -R $mm10min    -o /tmp/$child.dam.allele_depth    -I /tmp/$dam.list --printBaseCounts    -L /tmp/$child.bed
	pid2=$!
	java -jar $GATK3446    -T DepthOfCoverage -R $mm10min    -o /tmp/$child.kid.allele_depth    -I /tmp/$kid.list --printBaseCounts    -L /tmp/$child.bed
	pid3=$!
	wait $pid1
	wait $pid2
	wait $pid3
done



java -jar /utils/appz/snpEff/SnpSift.jar filter "isVariant( GEN[$child] ) & isRef( GEN[$mother] ) & isRef( GEN[$father] ) & isHom( GEN[$mother] ) & isHom( GEN[$father] )"

java -jar /utils/appz/snpEff/SnpSift.jar filter "isVariant( GEN[$child] ) & isRef( GEN[$mother] ) & isRef( GEN[$father] ) & isHom( GEN[$mother] ) & isHom( GEN[$father] )"

java -jar /utils/appz/snpEff/SnpSift.jar filter "isVariant( GEN[$child] ) & isRef( GEN[$mother] ) & isRef( GEN[$father] ) & isHom( GEN[$mother] ) & isHom( GEN[$father] )"

java -jar /utils/appz/snpEff/SnpSift.jar filter "isVariant( GEN[$child] ) & isRef( GEN[$mother] ) & isRef( GEN[$father] ) & isHom( GEN[$mother] ) & isHom( GEN[$father] )"



# IF KID IS HET && MOM & DAD ARE HOM REF

# IF KID IS HOM ALT && MOM & DAD ARE HOM REF

# IF KID IS HOM ALT && MOM IS HET && DAD IS HOM REF

# IF KID IS HOM ALT && MOM IS HOM REF && DAD IS HET

# IF KID IS HOM REF && MOM IS HOM ALT && DAD IS HOM REF

# IF KID IS HOM REF && MOM IS HOM REF && DAD IS HOM ALT










# FITLER OUT SEG DUPS (UCSC)
