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
    log: 
        out="../Logs/Analysis/stats.out.log",
        err="../Logs/Analysis/stats.err.log"
    conda: "../Envs/stats.yaml"
    shell:
        """
        echo 'Finished mapping' >> time.txt 
        date +%s%N >> time.txt 
        python Scripts/msGBS_STATS.py -i {input.bamOut} -o {output.statscsv} > {log.out} 2> {log.err}
        """

rule filter:
    input: 
        statscsv=expand("{path}/stats/stats.csv",path=config["output_dir"])
    log: 
        out="../Logs/Analysis/filter.out.log", 
        err="../Logs/Analysis/filter.err.log",
        log="../Logs/Analysis/filter.log.log"
    output:
        Data1=expand("{path}/stats/Data_1_Clusters_Target_vs_Reason_to_remove_8_15_1000_summed_per_species.txt", path=config["output_dir"]),
        Data2=expand("{path}/stats/Data_2_Clusters_filtered_due_to_homology_to_8_15_1000.txt", path=config["output_dir"]),
        Data3=expand("{path}/stats/Data_3_READ_COUNT_removed_CLUSTERS_SUM_8_15_1000.csv", path=config["output_dir"]),
        Data4=expand("{path}/stats/Data_4_SUM_8_15_1000.csv", path=config["output_dir"]),
        Data5=expand("{path}/stats/Data_5_SUM_MINREAD_FILTER_8_15_1000.csv", path=config["output_dir"])
    shell:
        """
        python Scripts/Parse_csv_Jelle_final.py -i {input.statscsv} -log {log.log} -f1 8 -f2 15 -f3 1000 -op ../Output/stats/ > {log.out} 2> {log.err}
        """

"""
rule logging:
    input:
        #alle logs
    output:
        #multiqc log
    conda:
        # multiqc log env
    shell:
        # multiqc command
        """