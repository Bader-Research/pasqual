# PASQUAL

Our assembler PASQUAL, short for PArallel SeQUence AssembLer, is
designed for shared memory parallelism, using OpenMP due to its good
tradeoff between performance and programmer productivity.


README
------

This file contains information about how to install and use PASQUAL, our parallel
de novo sequence assembler.


1) Install
----------

- Open the Makefile in a text editor. Adapt the settings to your environment, in
  particular choose your favorite compiler. Save and close.
  
- For the machines without SSE support, find "WITH_STNI = 1" in the Makefile
  and change it to "WITH_STNI = 0".

- Type make. This builds the object and binary files with default settings.

  Additional flags:
  - INDEX_LENGTH=<U32,U64> to specify the index data type size. Default: U32
                 When the sequence data sets are smaller than 4GB, use U32
                 to save memory usage. Otherwise, use U64.					
  - TIMING=<0,1> to specify if timing routines are enables or not. Default: 1


2) Run
------

Example workflow (adapt to your needs):
     
  iii) Assemble reads contained in reads.fa (read length 50) with a minimum overlap length of 31.
    
       $> bin/pasqual_omp_u32 or bin/pasqual_omp_u64 --singleread reads.fa -l 50 -t 31
       
       Try pasqual_omp_u32/pasqual_omp_u64 --help for more information.

   iv) Contigs are written to output file 'contigs.fa' in FASTA format.
