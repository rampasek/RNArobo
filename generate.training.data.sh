chromozones=`ls /projects/genomes/hg19/fa/`
chromozones="chrX.fa"

for des in `ls descriptors/`; do
	des2=`echo $des | sed s/.des$//g`
	for ch in $chromozones; do
		ch2=`echo $ch | sed s/.fa$//g`
		./rnarobo --tonly -q -c descriptors/$des /projects/genomes/hg19/fa/$ch > regression/$des2.$ch2.data.out
	done
done

