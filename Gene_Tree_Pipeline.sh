#!/bin/bash
#The job should run on the standard partition
#SBATCH -p standard
#The name of the job is tree
#SBATCH -J tree
#SBATCH --account=tyjames1
#The job requires 1 compute node
#SBATCH -N 1
#The job requires 1 task per node
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100gb
#The maximum walltime of the job is 14 days
#SBATCH -t 14-00:00:00
#SBATCH --mail-user=saleh_rahimlou@hotmail.com
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# Using a custom Python script to align DNA sequences with some assemblies
# Load the IQtree environment
# easel, mafft, and trimal packages are installed inside the IQtree environment
# You can change the gap treatment threshold and number of bootstraps
# Loading Python and hmmer modules are necessary

module load Bioinformatics
module load hmmer
module load python

./gene_extractor.py -seqfile=Saccharomyces.ERR1726578.augustus.abinitio.fasta -in_path=./Saccharomyces_assemblies -out_path=gene_alignment

# Removing empty files
for x in NODE_181_*.fasta;
do
        basename=$(echo $x | sed 's/'.fasta'//')
        cd $basename
        find . -type f -empty -print -delete
        cd ..
done

# Converting aligned stk files to unaligned fasta files
for x in NODE_181_*/;
do
        cd $x
        for file in *.stk
        do
                basename=$(echo $file | sed 's/[\.].*//')
                echo $basename
                esl-reformat fasta "$file" > "$basename".fasta
        done
        cd ..
done

# Get only the first hit in the output files
for x in NODE_181_*.fasta;
do
        basename=$(echo $x | sed 's/'.fasta'//')
        cd $basename
        for i in *.fasta;
        do
                outname=$(echo $i | sed 's/'.fasta'//')
                outname+="_1.fasta"
                ~/bbmap/reformat.sh in=$i out=$outname reads=1
        done
        cd ..
done

# Change the header of the sequence to the genome assembly accession
for x in NODE_181_*.fasta;
do
        basename=$(echo $x | sed 's/'.fasta'//')
        cd $basename
        for i in *_1.fasta;
        do
                header=$(echo $i | sed 's/'_1.fasta'//')
                filename=$(echo $i | sed 's/'_1.fasta'//')
                cat $i | sed 's/^>.*$/>'$header'/' > "$filename"_2.fasta
        done
        cd ..
done

# Copy query fasta files to corresponding subject directory
for x in NODE_181_*.fasta;
do
        basename=$(echo $x | sed 's/'.fasta'//')
        cp $x $basename/
done

# Edit the query file name
for x in NODE_181_*.fasta;
do
        basename=$(echo $x | sed 's/'.fasta'//')
        cd $basename
        mv $x "$basename"_2.fasta
        cd ..
done

# Pool query and subject files in all directories
for x in NODE_181_*.fasta;
do
        basename=$(echo $x | sed 's/'.fasta'//')
        cd $basename
        basename+='_pooled.fasta'
        cat *_2.fasta > $basename
        cd ..
done

# Use MAFFT to align the pooled file & trim using Trimal and constructing the phylogenetic tree with IQtree
for x in NODE_181_*.fasta;
do
        basename=$(echo $x | sed 's/'.fasta'//')
        cd $basename
        align_input=$basename
        align_input+='_pooled.fasta'
        align_output=$basename
        align_output+='_mafft.fasta'
        trimal_output=$basename
        trimal_output+='_trimal.fasta'
        mafft --auto $align_input > $align_output
        trimal -in $align_output -out $trimal_output -gt 0.1
        iqtree -s $trimal_output --prefix $basename -m MFP -T AUTO
        cd ..
done

