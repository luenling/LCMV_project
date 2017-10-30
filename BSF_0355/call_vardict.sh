


$VARDICT/VarDict_java -G $REFGENOME -f 0.01 -b ../BQSR/BSF_0355_S10_real_viterbi_IDQS_bqsr.bam -c 1 -S 2 -E 3 -z 1 $BASEDIR/References/viruses_short.bed | $VARDICT/VarDict/teststrandbias.R | $VARDICT/VarDict/teststrandbias.R -N S10 -E -f 0.01 > test_S10.vcf
