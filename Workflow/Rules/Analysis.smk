# Authors v1.0 (Legacy):    ..., ..., ... & ...
# Authors v2.0 (New):       Jelle van der Heide
# Date:                     ../../..
# 
# The rules in this file handle analyzing mapper output, as well as creating a summary of log-files.
# They format mapper output in a way where a Python file can accept it, and summarize it into a csv-file
# this csv-file can in turn be parsed by another Python-file, which outputs another five analysis files.
# --------------------------------------------------------------------------------------------------------------------



# Before the bam-files can be analyzed properly, they are required to have readgroups in their headers,
# as the Python file is unable to parse them properly in it's current state. By using picard readgroups
# are added retroactively.
# -----     
# Input:    - Mapper output bam-files, lacking readgroups in their headers.
# Output:   - Mapper output bam-files, now with readgroups added to their headers.
rule bam_rg:
    params:
        sample='{sample}', 
        mapper=MAPPER
    input:
        bamIn=expand("{bam_dir}/{{mapper}}/mapping_sq_{{sample}}.bam",bam_dir=config["bam_dir"])
    output:
        bamOut=expand("{bam_dir}/{{mapper}}/mapping_rg_{{sample}}.bam",bam_dir=config["bam_dir"])
    log: 
        "../Logs/Analysis/bam_rg_{sample}_{mapper}.log"
    benchmark: 
        "../Benchmarks/bam_rg_{sample}_{mapper}.benchmark.tsv"
    conda:
        "../Envs/bam_rg.yaml"
    #threads: NULL
    shell:
        """
        picard AddOrReplaceReadGroups \
            I={input.bamIn} \
            O={output.bamOut} \
            RGLB={params.sample} \
            RGPL=illumina \
            RGPU=unit \
            RGSM={params.sample} \
            RGID={params.sample} \
            2> {log}
        """

# As the Python file accepts a single file to parse, and the mappers output multiple bam-files,
# the mapper output bam-files should be merged into a single file. By using samtools merge they
# are merged into a single file.
# -----     
# Input:    - Multiple readgrouped bam-files from one of the mappers.
# Output:   - A single readgrouped bam-file from one of the mappers.
rule bam_merge:
    params: 
        mapper=MAPPER
    input: 
        bamIn=expand("{bam_dir}/{{mapper}}/mapping_rg_{sample}.bam",bam_dir=config["bam_dir"], sample=SAMPLES)
    output:
        bamOut=expand("{analysis_out}/{{mapper}}/mapping_merged.bam",analysis_out=config["analysis_out"])
    log: 
        "../Logs/Analysis/bam_merge_{mapper}.log"
    benchmark: 
        "../Benchmarks/bam_merge_{mapper}.benchmark.tsv"
    conda: 
        "../Envs/bam_merge.yaml"
    #threads: NULL
    shell:
        """
        samtools merge - {input.bamIn} | 
        samtools sort > {output.bamOut} \
            2> {log}
        """

# This rule essentially just calls a Python file that is able to analyze a readgrouped bam-file.
# This file parses the bam-file, and moulds its contents into a human readable csv-file that
# summarizes how many reads were mapped to which read from the meta-reference from which sample. 
# -----     
# Input:    - A readgrouped bam-file containing all data originating from the bam-files from one of the mappers.
# Output:   - A human readable csv-file summarizing how many reads were mapped to which meta-reference
#             read from which sample.
rule stats:
    params:
        mapper=MAPPER
    input:
        bamIn=expand("{analysis_out}/{{mapper}}/mapping_merged.bam",analysis_out=config["analysis_out"])
    output:
        statscsv=expand("{analysis_out}/{{mapper}}/stats.csv",analysis_out=config["analysis_out"])
    log: 
        out="../Logs/Analysis/stats_{mapper}.out.log",
        err="../Logs/Analysis/stats_{mapper}.err.log"
    benchmark:
        "../Benchmarks/stats_{mapper}.benchmark.tsv"
    conda: 
        "../Envs/stats.yaml"
    #threads: NULL
    shell:
        """
        echo 'Finished mapping with {params.mapper}'  >> time.txt 
        date +%s%N >> time.txt 
        python Scripts/Stats.py \
            -i {input.bamIn} \
            -o {output.statscsv} \
            > {log.out} \
            2> {log.err}
        """

# This rule simply calls another Python file as well. This time around it analyzes the csv-file
# created by the previous rule, resulting in five files containing additional information for users
# to properly interpret their data and results by filtering the csv according to the filter configuration.
# -----     
# Input:    - The aforementioned human readable cssv-file summarizing how many reads were mapped 
#             to which meta-reference read from which sample. 
# Output:   - The files discussing or displaying specific data after applying user-defined filters to the csv file.
rule filter:
    params: 
        mapper=MAPPER,
        outprefix=expand("{analysis_out}/{{mapper}}/",analysis_out=config["analysis_out"])
    input: 
        statscsv=expand("{analysis_out}/{{mapper}}/stats.csv",analysis_out=config["analysis_out"])
    output:
        Data1=expand("{analysis_out}/{{mapper}}/Data_1_Clusters_Target_vs_Reason_to_remove_8_15_1000_summed_per_species.txt", analysis_out=config["analysis_out"]),
        Data2=expand("{analysis_out}/{{mapper}}/Data_2_Clusters_filtered_due_to_homology_to_8_15_1000.txt", analysis_out=config["analysis_out"]),
        Data3=expand("{analysis_out}/{{mapper}}/Data_3_READ_COUNT_removed_CLUSTERS_SUM_8_15_1000.csv", analysis_out=config["analysis_out"]),
        Data4=expand("{analysis_out}/{{mapper}}/Data_4_SUM_8_15_1000.csv", analysis_out=config["analysis_out"]),
        Data5=expand("{analysis_out}/{{mapper}}/Data_5_SUM_MINREAD_FILTER_8_15_1000.csv", analysis_out=config["analysis_out"])
    log: 
        out="../Logs/Analysis/filter_{mapper}.out.log", 
        err="../Logs/Analysis/filter_{mapper}.err.log",
        log="../Logs/Analysis/filter_{mapper}.log.log"
    benchmark:
        "../Benchmarks/filter_{mapper}.benchmark.tsv"
    #conda: NULL
    #threads: NULL
    shell:
        """
        python Scripts/Parse_csv.py \
            -i {input.statscsv} \
            -log {log.log} \
            -f1 8 \
            -f2 15 \
            -f3 1000 \
            -op {params.outprefix} \
            > {log.out} \
            2> {log.err}
        """

# Finally, this rule attempts to summarize the results/execution success of a number of rules in this pipeline
# by using MultiQC over all generated logs. Granted, this skips over the majority of these logs, as 
# most of them are irrelevant, empty or perhaps only contain log-data originating from rules
# that generate logs that MultiQC is unable to parse.
# -----     
# Input:    - The files discussing or displaying specific data after applying user-defined filters to the csv file.
#             [NOTE] These files are not actually used by this rule, and are only defined as input because 
#                    the existence of these files confirms the completion of the entire pipeline, which would in turn
#                    mean all log-files are available to parse at that moment.
# Output:   - A MultiQC output file in HTML format, summarizing all data it could parse from the entirity
#             of the log folder at the time of executing the command. 
rule logging:
    params:
        logdir=config["log_dir"]
    input:
        Data1=expand("{analysis_out}/{mapper}/Data_1_Clusters_Target_vs_Reason_to_remove_8_15_1000_summed_per_species.txt", analysis_out=config["analysis_out"], mapper=MAPPER),
        Data2=expand("{analysis_out}/{mapper}/Data_2_Clusters_filtered_due_to_homology_to_8_15_1000.txt", analysis_out=config["analysis_out"], mapper=MAPPER),
        Data3=expand("{analysis_out}/{mapper}/Data_3_READ_COUNT_removed_CLUSTERS_SUM_8_15_1000.csv", analysis_out=config["analysis_out"], mapper=MAPPER),
        Data4=expand("{analysis_out}/{mapper}/Data_4_SUM_8_15_1000.csv", analysis_out=config["analysis_out"], mapper=MAPPER),
        Data5=expand("{analysis_out}/{mapper}/Data_5_SUM_MINREAD_FILTER_8_15_1000.csv", analysis_out=config["analysis_out"], mapper=MAPPER)
    output:
        alllogs=expand("{logsummary_dir}/multiQClogsummary.html", logsummary_dir=config["logsummary_dir"])
    log:
        out="../Logs/Analysis/logging.out.log",
        err="../Logs/Analysis/logging.err.log"
    benchmark:
        "../Benchmarks/logging.benchmark.tsv"
    conda:
        "../Envs/logging.yaml"
    #threads: NULL
    shell:
        """
        multiqc {params.logdir} \
            -n {output.alllogs} \
            > {log.out} \
            2> {log.err}
        """