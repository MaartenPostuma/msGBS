rule bam_rg:
    params:
        sample='{sample}'
    input:
        bamIn=expand("{path}/mapping/mapping_sq_{{sample}}.bam",path=config["output_dir"])
    output:
        bamOut=temp(expand("{path}/mapping/mapping_rg_{{sample}}.bam",path=config["output_dir"]))
    threads: 1
    conda:
        "../Envs/bam_rg.yaml"
    shell:
        """
        picard AddOrReplaceReadGroups I={input.bamIn} O={output.bamOut} RGLB={params.sample} RGPL=illumina RGPU=unit RGSM={params.sample} RGID={params.sample}
        """

rule bam_merge:
    input: 
        # bowtie2
        bamIn=expand("{path}/mapping/mapping_rg_{sample}.bam",path=config["output_dir"], sample=SAMPLES)
        # bwa-mem
        #bamIn=expand("{path}/mapping/mapping_rg_{sample}.bam",path=config["output_dir"], sample=SAMPLES)
    output:
        bamOut=temp(expand("{path}/mapping/mapping_merged.bam",path=config["output_dir"]))
    threads:1
    conda: "../Envs/bam_merge.yaml"
    shell:
        """
        date +%s%N >> time.txt
        samtools merge - {input.bamIn} | samtools sort > {output.bamOut}
        """

rule stats:
    input:
        bamOut=expand("{path}/mapping/mapping_merged.bam",path=config["output_dir"])
    output:
        statscsv=expand("{path}/stats/stats.csv",path=config["output_dir"])
    conda: "../Envs/stats.yaml"
    shell:
        "python Scripts/msGBS_STATS.py "
        "-i {input.bamOut} "
        "-o {output.statscsv}"
