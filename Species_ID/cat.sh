for f in *;
do
  	sed -e "s/\(^@M0.*\) .*$/\1;sample=${f%.*};/" $f \
        >> ../SL_mimam_ALL_CATout.fastq;
done
