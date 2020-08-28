# DEEPL
Introduction
============

DEEPL: 

    DEEP is a RNA-seq aligner for third generation sequencing (TGS) RNA-seq. The popurse of DEEP is to find out more information about transcripts from RNA-seq data. DEEP mapping fasta or fastq RNA-seq to reference with fasta formate. Reference files will contain a whole genome or more genomes, it also can be a transcipt or only a peiece of sequence with fasta formate. DEEP also accept other aligners primary align results as input for deep mode align. DEEP primary align step need a BWT index of reference, we will use BWA build BWT index. After that DEEP find overlapped max exact match (MEM) seeds of reads with BWT index, and find acceptable seed combinations, and achieve align results. DEEP output align result in SAM format, enabling interoperation with a large number of other tools that use SAM. At present, DEEP can only run on Linux operating system.

Version
============

    1.0.0


Memory usage
============
    
    DEEPL can align RNA-seq which sequencing from whole human transcripts to hg19 human reference within 16GB RAM memory. 

Installation
============

    cd DEEPL/src
    make

Commands and options
============

    Download BWA:
	 https://sourceforge.net/projects/bio-bwa/files/
    [build BWT index]:
        ./bwa index -p <bwa_index> ref.fa
        
    [build index]: 
        ./DEEPL index <input_fasta_file_route> <output_index_route>
    [align]: 
        ./DEEPL complete <-p/-f> -B <bwa_index> -H <output_index_route> -1 <input_file> -O <output_file> <other options>

    [options]:

        -B <STRING>:   reference BWT index path
        -H <STRING>:   hash files path, each need file's name must be *.hash with *.ann files
        -f         :   fa formate
        -q         :   fq formate
        -1 <STRING>:   input files, need *.fa,*.fq files,split with ','
        -O <STRING>:   output file,need enough space

        -r <INT>:      search area(500000)
        -p <INT>:      thread(1)
        -m <INT>:      match score(2)
        -s <INT>:      miss penalty score(4)
        -g <INT>:      gap penalty score(2)
        -d :           deep align mode
        -t <INT>:      cut tail exon shorter than <INT> (0)

Evaluation
============

     The data simulation and benchmarking scripts are available at: https://github.com/hitbc/deSALT.
     The source code of DEEPL are available at: https://github.com/cathy-houli/DEEPL
     evaluation program for real data are available at: https://github.com/cathy-houli/DEEPL/evaluation

Contact
============

    For advising, bug reporting and requiring help, please contact ye.wenzhu@gmail.com

