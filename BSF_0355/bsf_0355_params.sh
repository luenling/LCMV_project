#!/bin/bash
# author Lukas Endler
# created: 10/30/2017, 9:17:18 AM
# Changed: Time-stamp: <10/30/2017, 9:21:31 AM>
# Variables for scripts, looks for different users (vetgrid01 or vetlinux01)
shopt -s extglob

BASEDIR=$(dirname $BASH_SOURCE)
BASEDIR=${BASEDIR%\/LCMV_project\/BSF_0355}

RUN_ID=0355
SCRIPTS=${BASEDIR}/LCMV_project/BSF_0355/
REFGENOME=$BASEDIR/References/viruses_short.fasta
REFGENOME_FULL=$BASEDIR/References/viruses_Mus_musculus.GRCm38.fa.gz
PRIMERS=$BASEDIR/References/primers.fna

# setup everything for vetgrid01

if [ $USER == "vetgrid01" ]; then
  # general tools
  PICARD=/usr/local/Cellar/picard-tools/2.5.0/share/java/picard.jar
  SAMTOOLS=/usr/local/bin/samtools
  BWA=/usr/local//Cellar/bwa/0.7.15/bin/bwa
  BCFTOOLS=/usr/local/bin/bcftools
  # Lofreq executable
  LOFREQ=$BASEDIR/Tools/lofreq_star-2.1.2/bin/lofreq
  #varscan2
  VARSCAN=$BASEDIR/Tools/varscan/VarScan.v2.4.0.jar
  #VarDict java base dir
  VARDICT=$BASEDIR/Tools/VarDictJava/
  # setup variables for SNPEFF
  SNPEFF_DATA="/usr/local/Cellar/snpeff/4.2/share/snpeff/data"
  SNPEFF_CONF="/usr/local/Cellar/snpeff/4.2/share/snpeff/snpEff.config"
  SNPEFF="snpEff"
  SNPSIFT="SnpSift"
  # setup GATK path
  GATK=$BASEDIR/Tools/GenomeAnalysisTK_3.5.0.jar
  # setup path to BBDUK
  BBDUK=$BASEDIR/Tools/bbmap/bbduk.sh
  BBMERGE=$BASEDIR/Tools/bbmap/bbmerge.sh
  ADAPTS=$BASEDIR/Tools/bbmap/resources/adapters.fa
  # VPHASER2
  VPHASER2=$BASEDIR/Tools/viral-ngs/intrahost_alt2.py
  BAMTOOLS23=$BASEDIR/Tools/viral-ngs/tools/tools/binaries/V-Phaser-2.0/bamtools-2.3.0/lib

elif [ $USER == "vetlinux01" ] ; then
  # general tools
  PICARD=~/.linuxbrew/Cellar/picard-tools/2.12.1/share/java/picard.jar
  SAMTOOLS=samtools
  BWA=/home/vetlinux01/.linuxbrew/bin/bwa
  BCFTOOLS=/home/vetlinux01/bin/bcftools
  # Lofreq executable
  LOFREQ=$BASEDIR/Tools/lofreq_star-2.1.2/bin/lofreq
  #varscan2
  VARSCAN=$BASEDIR/Tools/VarScan.v2.4.2.jar
  #VarDict java base dir
  VARDICT=$BASEDIR/Tools/VarDictJava/
  # setup variables for SNPEFF
  SNPEFF_DATA="/usr/local/Cellar/snpeff/4.3p/share/snpeff/data"
  SNPEFF_CONF="/usr/local/Cellar/snpeff/4.3p/share/snpeff/snpEff.config"
  SNPEFF="snpEff"
  SNPSIFT="SnpSift"
  # setup GATK path
  GATK=~/Tools/GenomeAnalysisTK.jar
  # setup path to BBDUK
  BBDUK=~/Tools/bbmap_37.57/bbduk.sh
  BBMERGE=~/Tools/bbmap_37.57/bbmerge.sh
  ADAPTS=~/Tools/bbmap_37.57/resources/adapters.fa
  # VPHASER2
  VPHASER2=$BASEDIR/Tools/viral-ngs/intrahost_alt2.py
  BAMTOOLS23=$BASEDIR/Tools/viral-ngs/tools/tools/binaries/V-Phaser-2.0/bamtools-2.3.0/lib

fi


exit
