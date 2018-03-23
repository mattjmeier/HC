#!/bin/bash
# Takes gVCF outputs MVs

echo "Where is the file containing a list of offspring to test?"
read -e sample_list

echo "The following offspring will be used to generate raw Mendelian violation lists:"
for i in $(cat $sample_list)
do
echo $i
done

echo "The following directory will be used to search for merged gVCFs:"
rootdir="~/DNA/pedigree_BaP/tvc_output/2017_additional_offspring"
echo $rootdir
# read -e rootdir

read -r -p "Please ensure the above is correct before you continue!!! [y/N] " response
case $response in
    [yY][eE][sS]|[yY])

        echo "Going ahead with analysis..."
        # You have 1 second to kill it
        sleep 1

###########################################################
###########################################################
###  START COMMANDS #######################################
###########################################################
###########################################################

for kid in $(cat $sample_list)
do

echo "Generating MVs for $kid..."

# Get the merged gVCF
family_gVCF=$(find $rootdir -name "*$kid.trio.merged.sorted.vcf.gz" | grep -v archive)

echo "The location of the merged gVCF containing variant calls for each sample in the family:"
echo $family_gVCF

echo "Using GATK SelectVariants to output a comprehensive list of all Mendelian violations"
echo "Writing file to $rootdir/$kid/$kid.violations.raw.vcf"

# If family gVCF isn't indexed... index it
if [ ! -e $family_gVCF.tbi ]
then
	echo "Family gVCF not indexed; running tabix."
	tabix $family_gVCF
fi

screen -d -m -S "$kid".SelectVariants.MVs bash -c "
java -jar $GATK3446 -T SelectVariants -R $mm10min -V $family_gVCF -ped ~/DNA/pedigree_BaP/tvc_output/2017_additional_offspring/$kid/$kid.ped -mv  -o $rootdir/$kid/$kid.violations.raw.vcf

"

done
echo "Done."

###########################################################
###########################################################
###  END COMMANDS #########################################
###########################################################
###########################################################

        ;;
    *)
        echo "Quitting now"
        ;;
esac

