# MAG Reconstruction and annotation from experimental upwelling mesocosms
The steps and scripts that I used to analyze the MOSS Landing metagenomic data including both prokaryotic and eukaryotic lineages -- functions abundance and beyond

The data for this project were sequenced at the Princeton Sequencing core facility as 65bp PE fragments. A summary of the number of reads generated per sample are contained in this git "summary of reads per sample"

### 1. Demultiplex
    iu-demultiplex -s barcode_to_sample_8index.txt --r1 2101__Merge-sample-metagenomics-of-moss-landing-from-flowcell-HFCNGDMXY-on-2022-06-17_Read_1_passed_filter.fastq --r2 2101__Merge-sample-metagenomics-of-moss-landing-from-flowcell-HFCNGDMXY-
on-2022-06-17_Read_4_passed_filter.fastq --index 2101__Merge-sample-metagenomics-of-moss-landing-from-flowcell-HFCNGDMXY-on-2022-06-17_Read_2_Index_Read_passed_filter.fastq -o DEMULTIPLEX
#!/bin/bash

### 2. Set up the ini files for demultiplexing

    iu-gen-configs 00_DEMULTIPLEXING_REPORT
    ls *.ini | sed 's/\.ini//g' > x_samples.txt
    
### 3. Quality filter using the slurm scheduler

    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=1
    #SBATCH --time=02:00:00
    #SBATCH --mem=1Gb
    #SBATCH --array=1-15

    SAMPLE=$(sed -n "$SLURM_ARRAY_TASK_ID"p x_samples.txt)
    iu-filter-quality-minoche ${SAMPLE}.ini
    grep "@A00741" ${SAMPLE}-QUALITY_PASSED_R1.fastq | wc -l > ${SAMPLE}-QUALITY_PASSED_R1.fastq-reads.txt
    grep "@A00741" ${SAMPLE}-QUALITY_PASSED_R2.fastq | wc -l > ${SAMPLE}-QUALITY_PASSED_R2.fastq-reads.txt

### 4. Assemble the quality filtered reads.. There are a total of 3,287,724,146 read pairs used for the assembly.

    #!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=20
    #SBATCH --time=48:00:00
    #SBATCH --mem=500Gb

    megahit -1 B1D2T1_metagen-QUALITY_PASSED_R1.fastq,B1D3T1_metagen-QUALITY_PASSED_R1.fastq,B1D4T1A_metagen-QUALITY_PASSED_R1.fastq,B1D5T1A_metagen-QUALITY_PASSED_R1.fastq,B1D5T2A_metagen-QUALITY_PASSED_R1.fastq,B1D6T1B_metagen-QUALITY_PASSED_R1.fastq,B2D2T1_metagen-QUALITY_PASSED_R1.fastq,B2D3T1_metagen-QUALITY_PASSED_R1.fastq,B2D4T1A_metagen-QUALITY_PASSED_R1.fastq,B2D5T1A_metagen-QUALITY_PASSED_R1.fastq,B3D2T1_metagen-QUALITY_PASSED_R1.fastq,B3D3T1_metagen-QUALITY_PASSED_R1.fastq,B3D4T1A_metagen-QUALITY_PASSED_R1.fastq,B3D5T1A_metagen-QUALITY_PASSED_R1.fastq,B3D5T2C_metagen-QUALITY_PASSED_R1.fastq -2 B1D2T1_metagen-QUALITY_PASSED_R2.fastq,B1D3T1_metagen-QUALITY_PASSED_R2.fastq,B1D4T1A_metagen-QUALITY_PASSED_R2.fastq,B1D5T1A_metagen-QUALITY_PASSED_R2.fastq,B1D5T2A_metagen-QUALITY_PASSED_R2.fastq,B1D6T1B_metagen-QUALITY_PASSED_R2.fastq,B2D2T1_metagen-QUALITY_PASSED_R2.fastq,B2D3T1_metagen-QUALITY_PASSED_R2.fastq,B2D4T1A_metagen-QUALITY_PASSED_R2.fastq,B2D5T1A_metagen-QUALITY_PASSED_R2.fastq,B3D2T1_metagen-QUALITY_PASSED_R2.fastq,B3D3T1_metagen-QUALITY_PASSED_R2.fastq,B3D4T1A_metagen-QUALITY_PASSED_R2.fastq,B3D5T1A_metagen-QUALITY_PASSED_R2.fastq,B3D5T2C_metagen-QUALITY_PASSED_R2.fastq --k-list 39,59 -o x_FULL-megahit-assemblyk3959

### 5. Summarize the assembly using the quast
    
    #!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=1
    #SBATCH --time=00:25:00
    #SBATCH --mem=5Gb

    quast x_FULL-megahit-assembly-with-miseqreads/final.contigs.fa x_FULL-megahit-assemblyk3959/final.contigs.fa -o x_QUAST-assembly-comparison

### 6. Generate the anvio contigs database and annotate the DB with taxonomy and HMMs for virus etc.. 

    #!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=1
    #SBATCH --time=05:00:00
    #SBATCH --mem=10Gb

    #anvi-script-reformat-fasta /scratch/gpfs/WARD/JOE/MOSS_BLOOM/METAGENOMICS_20220617/DEMULTIPLEX/x_FULL-megahit-assemblyk3959/final.contigs.fa --simplify-names --min-len 2000 -o x_FULL-megahit-assemblyk3959/final.contigs.simplifiednames.fa --report-file x_FULL-megahit-assemblyk3959/final.contigs.simplifiednames.report.txt
    
    #anvi-gen-contigs-database -n FULL_megahit_assemblyk3959 -f /scratch/gpfs/WARD/JOE/MOSS_BLOOM/METAGENOMICS_20220617/DEMULTIPLEX/x_FULL-megahit-assemblyk3959/final.contigs.simplifiednames.fa -o /scratch/gpfs/WARD/JOE/MOSS_BLOOM/METAGENOMICS_20220617/DEMULTIPLEX/x_FULL-megahit-assemblyk3959/final.contigs.db -T 40
    
    #anvi-run-hmms -c /scratch/gpfs/WARD/JOE/MOSS_BLOOM/METAGENOMICS_20220617/DEMULTIPLEX/x_FULL-megahit-assemblyk3959/final.contigs.db -T 30

    #anvi-run-hmms -c /scratch/gpfs/WARD/JOE/MOSS_BLOOM/METAGENOMICS_20220617/DEMULTIPLEX/x_FULL-megahit-assemblyk3959/final.contigs.db -H /scratch/gpfs/WARD/JOE/ANVIO-HMM-DBs/HMMs_RNApol_A_and_B/HMM_RNA_b -T 30

    ## HMMS to detect giant viruses
    anvi-run-hmms -c /scratch/gpfs/WARD/JOE/MOSS_BLOOM/METAGENOMICS_20220617/DEMULTIPLEX/x_FULL-megahit-assemblyk3959/final.contigs.db -H /scratch/gpfs/WARD/JOE/EBAME-2023/03_Functions/HMM_Nitrogen_Fixation_TIGRFAM -T 30 --just-do-it
    anvi-run-hmms -c /scratch/gpfs/WARD/JOE/MOSS_BLOOM/METAGENOMICS_20220617/DEMULTIPLEX/x_FULL-megahit-assemblyk3959/final.contigs.db -H /scratch/gpfs/WARD/JOE/EBAME-2023/03_Functions/HMM_Photosynthesis -T 30 --just-do-it
    anvi-run-hmms -c /scratch/gpfs/WARD/JOE/MOSS_BLOOM/METAGENOMICS_20220617/DEMULTIPLEX/x_FULL-megahit-assemblyk3959/final.contigs.db -H /scratch/gpfs/WARD/JOE/EBAME-2023/01_Giant_Virus_markers/HMM_RNA_a -T 30 --just-do-it
    anvi-run-hmms -c /scratch/gpfs/WARD/JOE/MOSS_BLOOM/METAGENOMICS_20220617/DEMULTIPLEX/x_FULL-megahit-assemblyk3959/final.contigs.db -H /scratch/gpfs/WARD/JOE/EBAME-2023/01_Giant_Virus_markers/HMM_RNA_b-T 30 --just-do-it
    anvi-run-hmms -c /scratch/gpfs/WARD/JOE/MOSS_BLOOM/METAGENOMICS_20220617/DEMULTIPLEX/x_FULL-megahit-assemblyk3959/final.contigs.db -H /scratch/gpfs/WARD/JOE/EBAME-2023/01_Giant_Virus_markers/HMM_MCP -T 30 --just-do-it
    anvi-run-hmms -c /scratch/gpfs/WARD/JOE/MOSS_BLOOM/METAGENOMICS_20220617/DEMULTIPLEX/x_FULL-megahit-assemblyk3959/final.contigs.db -H /scratch/gpfs/WARD/JOE/EBAME-2023/01_Giant_Virus_markers/HMM_pATPase -T 30 --just-do-it
    anvi-run-hmms -c /scratch/gpfs/WARD/JOE/MOSS_BLOOM/METAGENOMICS_20220617/DEMULTIPLEX/x_FULL-megahit-assemblyk3959/final.contigs.db -H /scratch/gpfs/WARD/JOE/EBAME-2023/01_Giant_Virus_markers/HMM_Primase -T 30 --just-do-it
    anvi-run-hmms -c /scratch/gpfs/WARD/JOE/MOSS_BLOOM/METAGENOMICS_20220617/DEMULTIPLEX/x_FULL-megahit-assemblyk3959/final.contigs.db -H /scratch/gpfs/WARD/JOE/EBAME-2023/01_Giant_Virus_markers/HMM_TFIIS -T 30 --just-do-it
    anvi-run-hmms -c /scratch/gpfs/WARD/JOE/MOSS_BLOOM/METAGENOMICS_20220617/DEMULTIPLEX/x_FULL-megahit-assemblyk3959/final.contigs.db -H /scratch/gpfs/WARD/JOE/EBAME-2023/01_Giant_Virus_markers/HMM_VIRUS_149_hmms -T 30 --just-do-it
    anvi-run-hmms -c /scratch/gpfs/WARD/JOE/MOSS_BLOOM/METAGENOMICS_20220617/DEMULTIPLEX/x_FULL-megahit-assemblyk3959/final.contigs.db -H /scratch/gpfs/WARD/JOE/EBAME-2023/01_Giant_Virus_markers/HMM_VLTF3 -T 30 --just-do-it

    anvi-run-scg-taxonomy -c x_FULL-megahit-assemblyk3959/final.contigs.db -T 40

### Run the mapping of each sample back to the full assembly. Some of the parts of the steps below cannot be run at the same time. First you need to build the index for the assembly using bowtie2-build and this does not require an array, so you need to comment out the array line until this step finishes and then you can run the bowtie2 command. 

    #!/bin/bash
    #
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=10
    #SBATCH --time=04:00:00
    #SBATCH --mem=100Gb
    #SBATCH --array=1-15

    SAMPLE=$(sed -n "$SLURM_ARRAY_TASK_ID"p x_samples.txt)
    #bowtie2-build -f B1D4T1A-megahit/final.contigs.simplifiednames.fa B1D4T1A-megahit-contigs
    #bowtie2 -x B1D4T1A-megahit-contigs -1 ${SAMPLE}-QUALITY_PASSED_R1.fastq -2 ${SAMPLE}-QUALITY_PASSED_R2.fastq -S B1D4T1A-${SAMPLE}.sam
    #samtools view -bS -F 4 B1D4T1A-${SAMPLE}.sam > B1D4T1A-${SAMPLE}-RAW.bam
    #rm -rf B1D4T1A-${SAMPLE}.sam
    #anvi-init-bam B1D4T1A-${SAMPLE}-RAW.bam -o B1D4T1A-${SAMPLE}.bam
    #rm -rf B1D4T1A-${SAMPLE}-RAW.bam

### RUN CONCOCT automated binning using an estimted number of genomes at 20. Of course there are many more genomes than this, but we will generate MAGs by refining each of the metabins.

    #!/bin/bash

    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=20
    #SBATCH --time=10:00:00
    #SBATCH --mem=300Gb

    #cut_up_fasta.py /scratch/gpfs/WARD/JOE/MOSS_BLOOM/METAGENOMICS_20220617/DEMULTIPLEX/x_FULL-megahit-assemblyk3959/final.contigs.simplifiednames.fa -c 10000 -o 0 --merge_last -b x_full-contigs_10K.bed > x_full-contigs_10K.fa
    #concoct_coverage_table.py x_full-contigs_10K.bed FULL-megahit-assemblyk3959*.bam > x_full-coverage_table.tsv
    #concoct --composition_file x_full-contigs_10K.fa --coverage_file x_full-coverage_table.tsv -b x_full-concoct_output-20bins/ -t 40 -r 60 -c 20
    #merge_cutup_clustering.py x_full-concoct_output-20bins/clustering_gt1000.csv > x_full-concoct_output-20bins/clustering_merged.csv
    extract_fasta_bins.py /scratch/gpfs/WARD/JOE/MOSS_BLOOM/METAGENOMICS_20220617/DEMULTIPLEX/x_FULL-megahit-assemblyk3959/final.contigs.simplifiednames.fa x_full-concoct_output-20bins/clustering_merged.csv --output_path x_full-concoct_output-20bins/fasta_bins

#### Import the bins into the anvio contigs database. A file called CONCOCT-bins-for-anvio.txt is required for this step. I created this script using a bunch of manual unix babyscript. I know you can do this, even if takes a bit of time. The file just needs to have two columns and look like the example file in this git. "CONCOCT-bins-for-anvio.txt".

    #!/bin/bash

    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=20
    #SBATCH --time=00:10:00
    #SBATCH --mem=100Gb

    anvi-import-collection /scratch/gpfs/WARD/JOE/MOSS_BLOOM/METAGENOMICS_20220617/DEMULTIPLEX/x_full-concoct_output/CONCOCT-bins-for-anvio.txt -c x_FULL-megahit-assemblyk3959/final.contigs.db -p FULL_megahit_assemblyk3959-MERGE/PROFILE.db -C CONCOCT --contigs-mode

#### I ran the commands below directly on the head node, because I could not stand waiting for the slow scheduler! Each of these jobs took less than 1 min.
    anvi-get-sequences-for-hmm-hits -c x_FULL-megahit-assemblyk3959/final.contigs.db -p FULL_megahit_assemblyk3959-MERGE/PROFILE.db -B x_euk-bins-from-concoct20.txt --get-aa-sequences --hmm-sources HMM_RNA_a -C CONCOCT20 -o x_full-concoct_output-20bins/x_euk-bins-from-concoct20-HMM_RNAa.faa
    anvi-estimate-genome-completeness -c x_FULL-megahit-assemblyk3959/final.contigs.db -p FULL_megahit_assemblyk3959-MERGE/PROFILE.db -C CONCOCT20 -o x_full-concoct_output-20bins/x_concoct20-completion-estimates.txt
    anvi-estimate-scg-taxonomy -c x_FULL-megahit-assemblyk3959/final.contigs.db -p FULL_megahit_assemblyk3959-MERGE/PROFILE.db -C CONCOCT20 -o x_full-concoct_output-20bins/x_concoct20-taxonomy-estimates.txt

### Now you can summarize the CONCOCT Bins and decide which ones you want to refine and go from there using a manual approach to generate a clean MAG collection.

    anvi-summarize -c x_FULL-megahit-assemblyk3959/final.contigs.db -p FULL_megahit_assemblyk3959-MERGE/PROFILE.db -C CONCOCT20 -o FULL_megahit_assemblyk3959-MERGE-CONCOCT20-SUMMARY

#### Fire up a vis machine and run anvi refine on all of your metabins. This is how to get the refinement session started. Once you have finished binning, you can run the sammarize script again. 

    anvi-refine -c x_FULL-megahit-assemblyk3959/final.contigs.db -p FULL_megahit_assemblyk3959-MERGE/PROFILE.db -C CONCOCT20 -b bin_1
