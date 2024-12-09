rule bam_rg:
    params:
        sample='{sample}'
    input:
        bamIn=expand("{path}/mapping/mapping_sq_{{sample}}.bam",path=config["output_dir"])
    output:
        bamOut=expand("{path}/mapping/mapping_rg_{{sample}}.bam",path=config["output_dir"])
    threads: 1
    log: "../Logs/Analysis/bam_rg_{sample}.log"
    conda:
        "../Envs/bam_rg.yaml"
    shell:
        """
        picard AddOrReplaceReadGroups I={input.bamIn} O={output.bamOut} RGLB={params.sample} RGPL=illumina RGPU=unit RGSM={params.sample} RGID={params.sample} 2> {log}
        """

rule bam_merge:
    input: 
        bamIn=expand("{path}/mapping/mapping_rg_{sample}.bam",path=config["output_dir"], sample=SAMPLES)
    output:
        bamOut=expand("{path}/mapping/mapping_merged.bam",path=config["output_dir"])
    threads:1
    conda: "../Envs/bam_merge.yaml"
    log: "../Logs/Analysis/bam_merge.log"
    shell:
        """
        samtools merge - {input.bamIn} | samtools sort > {output.bamOut} 2> {log}
        """

rule stats:
    input:
        bamOut=expand("{path}/mapping/mapping_merged.bam",path=config["output_dir"])
    output:
        statscsv=expand("{path}/stats/stats.csv",path=config["output_dir"])
    conda: "../Envs/stats.yaml"
    shell:
        """
        echo 'Finished mapping' >> time.txt 
        date +%s%N >> time.txt 
        python Scripts/msGBS_STATS.py -i {input.bamOut} -o {output.statscsv}
        """
 # this also has a bit of outpout that could be sent to a log
