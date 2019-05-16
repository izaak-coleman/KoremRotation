# Builds a mantis data structure from input:
# 1: mantis exe path
# 2: squeakr exe path
# 3: kmer length
# 4: s parameter
# 5: output squeakr path
# 6: input list of datasets. The filename of the list of datasets should be the study title

mantis=$1
squeakr=$2
k=$3
s=$4
sqkr_out=$5
fq_list=$6
study_title=${fq_list::(-4)}

# Build squeakrs and .lst file (list of paths to squeakrs) from input fastqs
# If old .lst exists, touch and delete
touch $sqkr_out/${study_title}_k${k}s${s}.lst
rm $sqkr_out/${study_title}_k${k}s${s}.lst
for file in `cat $fq_list`; do
  fname=`basename $file`
  $squeakr count -e -k $k -s $s -t 1 -o $sqkr_out/${fname::(-6)}_k${k}s${s}.squeakr $file
  echo $sqkr_out/${fname::(-6)}_k${k}s${s}.squeakr >> $sqkr_out/${study_title}_k${k}s${s}.lst 
done

# build mantis
mkdir mantis_cqf_${study_title}_k${k}s${s}
$mantis build -s $s -i $sqkr_out/${study_title}_k${k}s${s}.lst -o mantis_cqf_${study_title}_k${k}s${s}

# Write metadata report
echo "K: $k" >> mantis_cqf_${study_title}_k${k}s${s}/build_report.txt
echo "S: $s" >> mantis_cqf_${study_title}_k${k}s${s}/build_report.txt
echo "fastq datasets:" >> mantis_cqf_${study_title}_k${k}s${s}/build_report.txt
cat $fq_list  >> mantis_cqf_${study_title}_k${k}s${s}/build_report.txt
echo "squeakrs built:" >> mantis_cqf_${study_title}_k${k}s${s}/build_report.txt
cat $sqkr_out/${study_title}_k${k}s${s}.lst >> mantis_cqf_${study_title}_k${k}s${s}/build_report.txt

