# VCF-Sample-Filter
A tool to efficiently subset very large multi-sample VCFs for specific samples. 
________________________________________________________________________________

Use at your own risk. I cannot provide support. All information obtained/inferred with this script is without any implied warranty of fitness for any purpose or use whatsoever.

ABOUT:

This program takes as input a compressed or uncompressed multi-sample VCF and also a list of samples. The program filters the VCF for the samples supplied nad generates a new VCF. This program was tested on a Linux workstation. 

RATIONAL: 

There are existing tools to filter VCFs, including bcftools. However, I ran into memory issues when filtering very large (multi-Gb) VCFs. This program aims to address this by using multi-threading and adjusting memory usage.  

PREREQUISITES:

1) C++ compiler (for example GCC/G++) 

INSTALLATION:

This program should work on all systems. Download the .cpp file and compile according to your system. For Linux Ubuntu, compile with g++ (g++ -std=c++11 -O3 -o VSF VCF_SampleFilter_V1_1.cpp -lz -lpthread).

COMMAND LINE OPTIONS: 

<b>-i</b> Input VCF (either compressed or uncompressed)

<b>-o</b> Output file name 

<b>-s</b> One column list of samples (must match exactly the sample names in the input VCF)

<b>-t</b> Number of threads (start with 2 for the initial test) 

<i>EXAMPLE USAGE</i>

./VSF -i input.vcf.gz -o output.vcf -s samples.txt -t 8
