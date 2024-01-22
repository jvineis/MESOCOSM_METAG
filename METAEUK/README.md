## These are the steps that I used to merge the mapping results to the taxonomy output from metaeuk. I applied this approach to each independent assembly and the co-assembly of all reads. I'm running what this looks like for the independent assemblies. However, this approach would look the same if you were to run it for a single assembly except you would not need the array related parts of the sbatch command.

### An example of an individual assembly of nearly 2 million PE reads using megahit looked like this.

    #!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=20
    #SBATCH --time=24:00:00
    #SBATCH --mem=500Gb
    megahit -1 B3D5T1A_metagen-QUALITY_PASSED_R1.fastq -2 B3D5T1A_metagen-QUALITY_PASSED_R2.fastq -o B3D5T1A_metagen-megahit

### Build an anvio contigs database for the assembly and export the prodigal gene calls.

    #!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=1
    #SBATCH --time=02:00:00
    #SBATCH --mem=5Gb
    anvi-script-reformat-fasta B1D6T1B_metagen-megahit/final.contigs.fa --simplify-names --min-len 2000 -o B1D6T1B_metagen_simplified_names.fa --report-file B1D6T1B_metagensimplifiednames.report.txt
    anvi-gen-contigs-database -f B1D6T1B_metagen_simplified_names.fa -o B1D6T1B_metagen-contigs.db
    anvi-export-gene-calls -c B1D6T1B_metagen-contigs.db -o B1D6T1B_metagen-prodigal-gene-calls.txt --gene-caller prodigal
    
### Metaeuk call all the genes and taxonomy for each of the assemlbies. The fasta file input to the easy-predict command is the simplified names file create above 

    #!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=1
    #SBATCH --time=08:00:00
    #SBATCH --mem=150Gb
    #SBATCH --array=1-15

    SAMPLE=$(sed -n "$SLURM_ARRAY_TASK_ID"p x_samples.txt)

    metaeuk easy-predict ${SAMPLE}_simplified_names.fa /projects/WARD/JOE/DBs/metaeuk/MMETSP_uniclust50_MERC ${SAMPLE}_predsResults ${SAMPLE}_tempFolder
    metaeuk taxtocontig ${SAMPLE}_tempFolder/latest/contigs ${SAMPLE}_predsResults.fas ${SAMPLE}_predsResults.headersMap.tsv /projects/WARD/JOE/DBs/metaeuk/MMETSP_zenodo_3247846_uniclust90_2018_08_seed_valid_taxids ${SAMPLE}_taxResults ${SAMPLE}_tmpDir --majority 0.5 --tax-lineage 1 --lca-mode 2

### Now you can create a file that contains the metaeuk annotation into anvio using the command below. This command will output a file for the metaeuk gene calls that can be concatenated with the 

    python ~/scripts/convert-metaeuk-to-anvio-genecalls.py --gff B2D4T1A_metagen_predsResults.gff --faa B2D4T1A_metagen_predsResults.fas --out B2D4T1A_metagen-metaeuk-external-gene-calls.txt
    
### Creeate the 
