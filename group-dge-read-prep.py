#!/usr/bin/python
#
# Connor Burbridge
# Global Institute for Food Security
# connor.burbridge@gifs.ca
#
##########################################################################################
# Setup
# script to help automate the dge-read-prep step on a grouped basis
# looks for a samples_config.tsv file formatted as:
# treatment1	read-file0-R1.fq	read-file0-R2.fq	read-file1-R1.fq...
# treatment2	read-file0-R1.fq	read-file2-R2.fq ...
# and so on
# two other files are also required for the script to run:
# one file that contains the BASE NAME of the control read files mentioned above
# i.e.
# read-file0-
# read-file1-
# and the same type of file for the treatment condition. Please make sure that the names
# of these files match the treatment name in the samples_config.tsv file.
# You should also ensure that you have Trimmomatic, FastQC, Hisat2, StringTie, R and the 
# DESeq2 R library installed and in your current path for things to run smoothly.
#
##########################################################################################
# Genome preparation
# Prior to running analysis you should have the genome indexed. You should also ahve the 
# appropriate gff file for the genome of interest. When you have these, you must update 
# the genome and gff paths in teh dge-read-preparation.sh script to the paths of your 
# files. The genome must me indexed using the program Hisat2.
#
##########################################################################################
# Usage
# To run this script, enter the directory that contains the samples_config.tsv, 
# treatment1.txt and treatment2.txt files and run the command below:
# Command: python group-dge-read-prep.py path-to-reads-dir threads treatment1-count treatment2-count
# path-to-dir: absolute path to the directory that contains all read files
# threads: number of threads to use during analysis
# treatment1-count: number of samples associated with treatment1
# treatment2-count: number of samples associated with treatment2
#
##########################################################################################
# This script performs some checks that prevents re-doing analysis when results files are
# are already present. Should the pipeline fail at any step, please remove the output files
# from the most recent step and run the same command again. The pipeline should pick up 
# at that step.
# Also, should you want to remove ALL directories and files generated by the pipeline,
# run the following command from the directory containing the samples_config.tsv and 
# treatment files.
# command: python group-dge-read-prep.py clean
# This command will PERMANENTLY delete all results generated from any previous analyses.
# Be careful!
# It is also recommended that you re-direct the STDOUT from this pipeline to a log file
# with a simple > log.txt (or >> log.txt to append) at the end of any command.
#
##########################################################################################


from subprocess import call
from subprocess import check_call
from subprocess import run
import fileinput, re, sys, shutil, os, argparse

###########################################################################
# callR(treatmentList, sampleCount1, sampleCount2)
# The final method in the differential gene expression pipeline. Calls
# the Rscript DESeq2.txt after generating transcript and gene count files.
# 3 arguments:
# treatmentList: python list of the treatments previosuly generated
# sampleCount1: the number of replicates in treatment 1
# sampleCount2: the number of replicates in treatment 2
###########################################################################
def callR(treatmentList, sampleCount1, sampleCount2):
        if "gene_count_matrix.csv" in os.listdir("./"):
                print("Skipping prepDE.py script. Count matrix already present in current directory.")
        else:
                call(["prepDE.py", "-i", "gtffiles.txt"])
        
        try:
                call(["mkdir", "deseq2-output"])
        except OSError as e:
                print("Error: %s - %s." % (e.filename, e.strerror))

        try:
                dgeLogFile = open("deseq2-log.txt", "w")
                call(["cd", "deseq2-output"])
                call(["/usr/bin/Rscript", "--vanilla" ,"/u1/cbe453/gene-expression/GIFS-DESeq2/DESeq2.txt", treatmentList[0], treatmentList[1], str(sampleCount1), str(sampleCount2), "./gene_count_matrix.csv"], stdout=dgeLogFile, stderr=dgeLogFile)
                dgeLogFile.close()
        except Exception as e:
                print("Unable to change directory and run DESeq2.txt...")
        sys.exit(0)

############################################################################
# prepGTF(treatmentList, outfile)
# This method generates gtffiles.txt file which contains paths to each
# gtf file generated by the dge-read-preparation.sh script.
# 2 arguments:
# treatmentList: python list of treatments previously generated
# outfile: an open file handle with hardcoded name "gtffiles.txt"
############################################################################
def prepGTF(treatmentList, outfile):
        gtffiles= []
        # Generate the file that contains the paths to the generated GTF files
        for treatment in treatmentList:
                treatmentPath = ("./" + treatment + "/stringtie-output/")	
                gtffiles = os.listdir(treatmentPath)
                cwd = os.getcwd()
                i = 0
                print(treatment)
                for gtf in gtffiles:
                        print(gtf)
                        i += 1
                        outfile.write(treatment + str(i) + "\t" + cwd + "/" + treatmentPath + gtf + "\n")

############################################################################
# clean()
# Simple method to remove any and all output generated by this pipeline.
# Will attempt to remove everythingg, even if it doesn't exist.
############################################################################
def clean():
        try:
                print("Flag --clean supplied, removing all output files and directories...")
                call(["rm", "gtffiles.txt"])
                call(["rm", "gene_count_matrix.csv"])
                call(["rm", "transcript_count_matrix.csv"])
                call(["rm", "-rf", "deseq2-output"])
       	        call(["rm", "Rplots.pdf"])
       	except Exception as e:
                print("Error: %s - %s." % (e.filename, e.strerror))
        sys.exit()

############################################################################
# dgeReadPrep(args, treatmentFile)
# This method calls dge-read-preparation.sh, which is a shell script that
# does most of the heavy lifting for this pipeline. It performs all steps
# up to and including the generation of the indiviual gtf files for each
# replicate.
# 2 arguments:
# args: the argparse representation of command line arguments
# treatmentList: file handle for the samples_config.tsv file
############################################################################
def dgeReadPrep(args, treatmentFile):
	# iterate through each treatment supplied in the samples_config.tsv file
        treatmentList = []
        for line in treatmentFile:
                splitLine = line.split()
                treatment = splitLine[0]
                treatmentList.append(splitLine[0])
                prepLogFile = open("./read-prep-log.txt", "w")

                if args.clean == "true":
                        try:
                                shutil.rmtree(treatment)
                        except Exception as e:
                                print("Error: %s - %s." % (e.filename, e.strerror))
                else:
                        print("Preparing reads for genotype: " + splitLine[0])
                        call(["mkdir", splitLine[0]])
                        check_call(["bash", "-i", "./dge-read-preparation.sh", treatment, str(args.threads), args.reads, args.genome, args.gff], stdout=prepLogFile, stderr=prepLogFile)

        if (args.clean == 'true'):
                clean()

        prepLogFile.close()
        treatmentFile.close()
        return treatmentList
        

#############################################################################
# main(args)
# The main method of the pipeline. Does some variable handling and calls all
# major sub-methods.
# 1 argument:
# args: argparse representation of command line arugments
#############################################################################
def main(args):

        # initialize conda env
        call(["conda", "init", "bash"])

        # check for the clean option and manage command-line aruguments
        if (args.clean == "true"):
                print("Option clean found. Removing generated directories for each treatment.")
        else:
                reads_dir = args.reads
                threads = args.threads
                sampleCount1 = args.countOne
                sampleCount2 = args.countTwo
        
        # IO mangement
        treatmentFile = open("samples_config.tsv", "r")
        treatmentList = []

        # Call the dgeReadPrep functino to perform steps up to and including
        # transcript counting
        treatmentList = dgeReadPrep(args, treatmentFile)

        # Call the prepGTF method to prepare gtffiles.txt
        outfile = open("gtffiles.txt", "w")
        prepGTF(treatmentList, outfile)
        outfile.close()

        # Call the final DESeq2.txt Rscript for differential gene expression calculation
        callR(treatmentList, sampleCount1, sampleCount2)
        return(1)

#############################################################################
# Start up method
# handles some simple argparse arguments and hands them off to the main
# method.
#############################################################################
if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='Main script to use when running differential gene expression pipeline. Make sure to downlaod the latest set of scripts from https://github.com/cbe453/GIFS-DESeq2.')
        parser.add_argument('--reads', type=str, dest='reads', help='full path to to the directory containing the desired read files.')
        parser.add_argument('--threads', type=int, dest='threads', help='number of desired threads for processes that support multi-threaded computation.')
        parser.add_argument('--sample_count_1', type=int, dest='countOne', help='number of replicates for the first treatment in the samples_config.tsv file.')
        parser.add_argument('--sample_count_2', type=int, dest='countTwo', help='number of replicates for the second treatment in the samples_config.tsv file.')
        parser.add_argument('--clean', type=str, dest='clean', help='if the clean flag is true, remove all output directories and files. Acceptable values are true or false. Default value is false.', default='false')
        parser.add_argument('--genome', type=str, dest='genome', help='path and prefex of the Hisat2 indices for the genome of interest.')
        parser.add_argument('--gff', type=str, dest='gff', help='path and name of the gff file for the genome of interest.')
        args = parser.parse_args()
        sys.exit(main(args))
