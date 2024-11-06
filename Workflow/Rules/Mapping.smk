#rule mapping_Star:
#    params:
#        tmpdirthis=tmpdirthis,
#        inputdir=expand("{path}/output_denovo",  path=config["output_dir"]),
#        outputdir=expand("{path}/mapping",  path=config["output_dir"]),
#        barcodes=expand("{path}/{bar}", path=config["input_dir"], bar=config["barcode_filename"]),
#    input:
#        refBlasted=expand("{path}/output_blast/Eukaryota_ref.fa",path=config["output_dir"]),
#        joined_All=expand("{path}/output_denovo/all.joined.fastq.gz",  path=config["output_dir"]),
#        merged_All=expand("{path}/output_denovo/all.merged.fastq.gz",  path=config["output_dir"])
#    output:
#        log=expand("{path}/mapping/mapping_variantcalling.log",path=config["output_dir"]),
#        bamOut=expand("{path}/mapping/out.bam",path=config["output_dir"])
#    threads: 40
#    shell:
#        "python Scripts/Map_STAR_snake.py --tmpdir {params.tmpdirthis} "
#        "--input_dir {params.inputdir} "
#        "--output_dir {params.outputdir} "
#        "--threads {threads} "
#        "--barcodes {params.barcodes}"


# --un <path> stuurt reads die niet alignen ergens heen
# --un-conc <path> doet dit voor read paren die niet aans;uitend mappen
# --met 30
rule mapping_Bowtie2: #denk hier aan de -k 10 optie. dit zou volgens mij maximaal 10 multimaps moeten geven    
    input:
        refBlasted=expand("{path}/output_blast/Eukaryota_ref.fa",path=config["output_dir"]),
        r1=expand("{path}/Preprocessing/demultiplexed_R1.fq.gz",  path=config["output_dir"]),
        r2=expand("{path}/Preprocessing/demultiplexed_R2.fq.gz",  path=config["output_dir"])
    conda: "../Envs/bowtie2.yaml"
    output:
        bamOut=expand("{path}/mapping/mapping.bam",path=config["output_dir"])
    shell:
        """
        bowtie2-build -f {input.refBlasted} ../Output/mapping/index
        bowtie2 -x ../Output/mapping/index -1 {input.r1} -2 {input.r2} -q --end-to-end --very-fast -k 10 --threads 10 -S ../Output/mapping/mapping.sam
        samtools view -b -o ../Output/mapping/mapping.bam ../Output/mapping/mapping.sam 
        """

rule mapping_BWA:

rule mapping_GEM3:

rule mapping_Segemehl:

