# Builds a mantis data structure from input:
# 1: mantis exe path
# 2: squeakr exe path
# 3: kmer length
# 4: s parameter
# 5: output squeakr path
# 7: input list of datasets. The filename of the list of datasets should be the study title

$mantis=$1
$squeakr=$2
$k=$3
$s=$4
$sqkr_out=$5
$mtis_out=$6
$fq_list=$7


# Build squeakrs and .lst file (list of paths to squeakrs) from input fastqs
for file in `cat $fq_list`; do
  $squeakr count -e -k $k -s $s -t 1 -o $sqkr_out/${file::(-5)}.squeakr $file
  echo $sqkr_out/${file::(-5)}.squeakr >> $sqkr_out/${fq_list::(-3)}.lst 
done

# build mantis
$mantis build -s $s -i $sqkr_out/${fq_list::(-3)}.lst -o $mtis_out
