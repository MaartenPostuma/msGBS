# Authors v1.0 (Legacy):    ..., ..., ... & ...
# Authors v2.0 (New):       Jelle van der Heide
# Date:                     ../../..
# 
# This file contains all msGBS pipeline functionality related to pre-processing raw reads in preparation for analysis. 
# Additionally, these rules can process a single pair of read files., as well as multiple, granted the barcode-file adheres
# to the guidelines stated in the README.
# --------------------------------------------------------------------------------------------------------------------



# All reads are trimmed for adapters and Poly-G's. Besides trimming, PCR duplicates and low-quality reads
# are filtered out  with fastP.
# -----
# Input:    - The read files with raw reads.
# Output:   - The files with trimmed reads, from which low quality reads and PCR duplicates have been removed. 
#           - A report in HTML format containing quality control results.
#           - A report in JSON format containing quality control results.
rule deduplicate_trim:
    params:
        run="{run}",
        adapter1=config["adapter1"],
        adapter2=config["adapter2"]
    input:
        reads1=expand("{input}/{{run}}_R1.fq.gz",input=config["input_dir"]),
        reads2=expand("{input}/{{run}}_R2.fq.gz",input=config["input_dir"])
    output:
        filtered1=temp(expand("{tmp_dir}/Preprocessing/Deduplicated/{{run}}_R1.fq.gz", tmp_dir=config["tmp_dir"])),
        filtered2=temp(expand("{tmp_dir}/Preprocessing/Deduplicated/{{run}}_R2.fq.gz", tmp_dir=config["tmp_dir"])),
        fastp_html=expand("{output_dir}/Preprocessing/Fastpreports/{{run}}Preprocessing.html", output_dir=config["output_dir"]),
        fastp_json=expand("{output_dir}/Preprocessing/Fastpreports/{{run}}Preprocessing.json", output_dir=config["output_dir"])
    log: 
        expand("{output_dir}/Logs/Preprocessing/deduplicate_trim_{{run}}.log", output_dir=config["output_dir"])
    benchmark: 
       "../Benchmarks/deduplicate_trim.benchmark_{run}.tsv"
    conda: 
        "../Envs/deduplication.yaml"
    resources:
        mem_mb= 10000,
        runtime= 120,
        cpus_per_task= 8
    threads: 
        8
    shell:
        """
        fastp --in1 {input.reads1} \
            --in2 {input.reads2} \
            --out1 {output.filtered1} \
            --out2 {output.filtered2} \
            --adapter_sequence {params.adapter1} \
            --adapter_sequence_r2 {params.adapter2} \
            --dedup \
            --trim_poly_g \
            --umi \
            --umi_loc per_read \
            --umi_len 3 \
            -j {output.fastp_json} \
            -h {output.fastp_html} \
            -w {threads} \
            2> {log}
        """ 

# Considering the barcode-file contains numerous columns of data that is irrelevant in this pipeline, this rule handles extracting 
# the relevant data from the file. The relevant data is stored in a new, split-up barcodefile which
# can in turn be used as input for other Rules instead of the bloated original barcodefile. Afterwards, stacks is used to process radtags, which can be fetched 
# from the filtered barcode-file. 
# [NOTE] As snakemake is somewhat strict in declaring which rules expext what data, this rule's output is officially considered to be a directory. However, 
#        The proper output should be considered to be the reads, demultiplexed in separate files named after the correct sample. As multiple runs can be processed
#        at once, they will be temporarily stored separately, to be merged or moved by the rule that follows this one. The temproary directory that this rule
#        outputs is not to be used, and is only to be considered a load bearing parameter for the actual output. 
# -----     
# Input:    - A .tsv-file containing read barcodes, which includes data surrounding their origin. 
#             [NOTE]The contents of this file should always follow the guidelines
#             stated in the README. The related parameters in the config-file should be adjusted accordingly.
#           - The quality filtered read files.
# Output:   - A filtered version of the before mentioned file. Data irrelevant to the pre-processing-steps
#             was removed from this version of the file to optimize data-extraction and reproducability.
#           - The demultiplexed read files 
#             [NOTE]As stated above, these files are not directly mentioned as output in the rule definition.
rule split_barcodes:
    input:
        barcodefile=expand("{input}/{barcodes}", input=config["input_dir"], barcodes=config["barcode_file"]),
    output:
        barcodefilefiltered=expand("{output_dir}/Preprocessing/Barcodesfiltered/{{run}}_barcodefiltered.tsv", output_dir=config["output_dir"]),
    log: 
        log1=expand("{output_dir}/Logs/Preprocessing/split_barcode_file_{{run}}.log", output_dir=config["output_dir"]),
    benchmark: 
       "../Benchmarks/split_barcode_file.benchmark_{run}.tsv"
    params:
        run="{run}"    
    resources:
        mem_mb= 1000,
        runtime= 2,
        cpus_per_task= 1
    conda:
        "../Envs/statsCombine.yaml"
    threads: 
        1
    shell:
         "Rscript Scripts/splitBarcode.R {input.barcodefile} {output.barcodefilefiltered} {params.run}"   



rule demultiplex:
    params:
        readdir=directory(expand("{tmp_dir}/Preprocessing/Demultiplexed", tmp_dir=config["tmp_dir"]))
    input:
        barcodefilefiltered=expand("{output_dir}/Preprocessing/Barcodesfiltered/{{run}}_barcodefiltered.tsv", output_dir=config["output_dir"]),
        filtered1=expand("{tmp_dir}/Preprocessing/Deduplicated/{{run}}_R1.fq.gz", tmp_dir=config["tmp_dir"]),
        filtered2=expand("{tmp_dir}/Preprocessing/Deduplicated/{{run}}_R2.fq.gz", tmp_dir=config["tmp_dir"])
    output:
        tmpdir=temp(directory(expand("{tmp_dir}/Preprocessing/Demultiplexed/{{run}}", tmp_dir=config["tmp_dir"])))
    log: 
        expand("{output_dir}/Logs/Preprocessing/process_radtags_{{run}}.log", output_dir=config["output_dir"])
    benchmark: 
       "../Benchmarks/split_barcode_file.benchmark_{run}.tsv"
    conda: 
        "../Envs/deduplication.yaml"
    resources:
        mem_mb= 10000,
        runtime= 240,
        cpus_per_task= 8
    threads: 
        8
    shell:    
        """
        mkdir {output.tmpdir}
        process_radtags \
            -1 {input.filtered1} \
            -2 {input.filtered2} \
            -b {input.barcodefilefiltered} \
            -o {output.tmpdir} \
            -r \
            --inline_inline \
            --renz_1 PacI \
            --renz_2 NsiI \
            --retain-header \
            --disable-rad-check \
            --threads {threads} 
            mv {output.tmpdir}/process_radtags*.log {log}
        """

# In case there number of samples exceeds the number of available barcodes, multiple runs may required, which results in multiple read-pair files
# However, as these runs re-use barcodes, process radtags is unable to properly process them all at once. This requires separate radtags results to be stored in separate 
# directories. This rule essentially undoes this separation by merging or moving sample files from each directory to one general directory. In case two
# files in different directories have a common name, they are to be considered the same sample and end up merged into a single file. Unique files are
# simply moved to the target directory.
# -----     
# Input:    - The demultiplexed read files. 
#             [NOTE]As stated above, these files are not directly mentioned as input in the rule definition. They are fetched from separate directories named after the run.
# Output:   - All demultiplexed read files, any duplicate sample merged into a single file. 
rule rename_samples:
    params:
        lambda w: {DEMULTIPLEXSAMPLES[w.demultiplexsamples]}, 
        readfile="{readfile}",
        tmp_dir=config["tmp_dir"],
        renamedunzipped=expand("{tmp_dir}/Preprocessing/samples/{{demultiplexsamples}}.{{readfile}}.fq", tmp_dir=config["tmp_dir"])
    input:
        tmpdir=expand("{tmp_dir}/Preprocessing/Demultiplexed/{run}", tmp_dir=config["tmp_dir"], run=RUN)
    output:
        renamed=expand("{output_dir}/Preprocessing/samples/nonmonos/{{demultiplexsamples}}.{{readfile}}.fq.gz",output_dir=config["output_dir"])
    log:
        expand("{output_dir}/Logs/Preprocessing/rename_samples_{{demultiplexsamples}}_{{readfile}}.log", output_dir=config["output_dir"])
    benchmark: 
       "../Benchmarks/rename_samples_{demultiplexsamples}_{readfile}.benchmark.tsv"
    #conda: NULL
    resources:
        mem_mb= 1000,
        runtime= 20,
        cpus_per_task= 1
    threads:
        1
    shell:
        """
        cat {params.tmp_dir}/Preprocessing/Demultiplexed/*/{wildcards.demultiplexsamples}.{params.readfile}.fq.gz > {output.renamed} 2> {log}
        """


rule move_monos_preprocess:
    input:
        monoReads_R1=expand("{reads_loc}/{{monos}}.1.fq.gz",reads_loc=config["reads_loc"]),
        monoReads_R2=expand("{reads_loc}/{{monos}}.2.fq.gz",reads_loc=config["reads_loc"])
    output:
        monoReads_R1=expand("{output_dir}/Preprocessing/samples/monos/{{monos}}.1.fq.gz",output_dir=config["output_dir"]),
        monoReads_R2=expand("{output_dir}/Preprocessing/samples/monos/{{monos}}.2.fq.gz",output_dir=config["output_dir"])
    threads: 
        1
    resources:
        mem_mb= 1000,
        runtime= 30,
        cpus_per_task= 1   
    shell:
        """
        cp {input.monoReads_R1} {output.monoReads_R1}
        cp {input.monoReads_R2} {output.monoReads_R2}
        """
rule move_nonmonos_preprocess:
    input:
        monoReads_R1=expand("{output_dir}/Preprocessing/samples/{{nonmonos}}.1.fq.gz",output_dir=config["output_dir"]),
        monoReads_R2=expand("{output_dir}/Preprocessing/samples/{{nonmonos}}.2.fq.gz",output_dir=config["output_dir"])
    output:
        monoReads_R1=expand("{output_dir}/Preprocessing/samples/nonmonos/{{nonmonos}}.1.fq.gz",output_dir=config["output_dir"]),
        monoReads_R2=expand("{output_dir}/Preprocessing/samples/nonmonos/{{nonmonos}}.2.fq.gz",output_dir=config["output_dir"]),
    threads: 
        1
    resources:
        mem_mb= 1000,
        runtime= 30,
        cpus_per_task= 1   
    shell:
        """
        cp {input.monoReads_R1} {output.monoReads_R1}
        cp {input.monoReads_R2} {output.monoReads_R2}
        """

rule createRef:
    input:
        indRef=expand("{ref_loc}/{monos}.fa",ref_loc=config["ref_loc"],monos=MONOS),
        monoReads_R1=expand("{output_dir}/Preprocessing/samples/monos/{monos}.1.fq.gz",output_dir=config["output_dir"],monos=MONOS),
        monoReads_R2=expand("{output_dir}/Preprocessing/samples/monos/{monos}.2.fq.gz",output_dir=config["output_dir"],monos=MONOS),
    output:
        ref=expand("{output_dir}/Blasting/Eukaryota_ref.fa",output_dir=config["output_dir"])
    threads: 
        1
    resources:
        mem_mb= 1000,
        runtime= 30,
        cpus_per_task= 1   
    shell:
        "cat {input.indRef} > {output.ref}"
        
rule mapping_Bowtie2_index:
    params: 
        indexprefix=expand("{output_dir}/Mapping/Index/Bowtie/index", output_dir=config["output_dir"])
    input: 
        refblasted=expand("{output_dir}/Blasting/Eukaryota_ref.fa",output_dir=config["output_dir"])
    output:
        index=expand("{output_dir}/Mapping/Index/Bowtie/index.1.bt2", output_dir=config["output_dir"])
    log:
        out=expand("{output_dir}/Logs/Mapping/mapping_bowtie2_index.out.log",output_dir=config["output_dir"]), 
        err=expand("{output_dir}/Logs/Mapping/mapping_bowtie2_index.err.log",output_dir=config["output_dir"])
    benchmark:
       "../Benchmarks/mapping_Bowtie2_index.benchmark.tsv"
    conda: 
        "../Envs/bowtie2.yaml"
    threads: 
        16
    resources:
        mem_mb= 10000,
        runtime= 480,
        cpus_per_task= 16 
    shell:
        """
        echo "Commencing Bowtie mapping" >> time.txt
        date +%s%N >> time.txt
        bowtie2-build \
            -f {input.refblasted} \
            {params.indexprefix} \
            -p {threads} \
            > {log.out} \
            2> {log.err}
        """



rule mapping_Bowtie:   
    params:
        sample='{sample}',
        multimap=config["multimap_bowtie"],
        indexprefix=expand("{output_dir}/Mapping/Index/Bowtie/index", output_dir=config["output_dir"])
    input:
        index=expand("{output_dir}/Mapping/Index/Bowtie/index.1.bt2", output_dir=config["output_dir"]),
        r1=expand("{output_dir}/Preprocessing/samples/{{type}}/{{sample}}.1.fq.gz",output_dir=config["output_dir"]),
        r2=expand("{output_dir}/Preprocessing/samples/{{type}}/{{sample}}.2.fq.gz",output_dir=config["output_dir"])
    output:
        samOut=temp(expand("{tmp_dir}/Mapping/Samout/Bowtie/{{type}}/mapping_sq_{{sample}}.sam",tmp_dir=config["tmp_dir"])),
        bamOut=temp(expand("{tmp_dir}/Mapping/Bamout/Bowtie/{{type}}/mapping_sq_{{sample}}.bam",tmp_dir=config["tmp_dir"]))
    log: 
        bowtie2=expand("{output_dir}/Logs/Mapping/{{type}}/mapping_bowtie2_{{sample}}_bt.log",output_dir=config["output_dir"]),
        samtools=expand("{output_dir}/Logs/Mapping/{{type}}/mapping_bowtie2_{{sample}}_st.log",output_dir=config["output_dir"])
    benchmark:
       "../Benchmarks/{type}/mapping_Bowtie2_{sample}.benchmark.tsv"
    conda: 
        "../Envs/bowtie2.yaml"
    threads: 
        8
    resources:
        mem_mb= 10000,
        runtime= 480,
        cpus_per_task= 8   
    shell:
        """
        bowtie2 \
            -x {params.indexprefix} \
            {params.multimap} \
            -1 {input.r1} \
            -2 {input.r2} \
            -q \
            --end-to-end \
            --very-fast \
            --threads {threads} \
            -S {output.samOut} \
            2> {log.bowtie2}
        samtools view \
            -b \
            -o {output.bamOut} {output.samOut} \
            2> {log.samtools}
        """

rule bam_rg_monos:
    params:
        sample='{sample}', 
        mapper=MAPPER
    input:
        bamIn=expand("{tmp_dir}/Mapping/Bamout/{{mapper}}/{{type}}/mapping_sq_{{sample}}.bam",tmp_dir=config["tmp_dir"])
    output:
        bamOut=expand("{output_dir}/Mapping/Bamout/{{mapper}}/{{type}}/mapping_rg_{{sample}}.bam",output_dir=config["output_dir"])
    log: 
        expand("{output_dir}/Logs/Analysis/{{type}}/bam_rg_{{sample}}_{{mapper}}.log",output_dir=config["output_dir"])
    benchmark: 
       "../Benchmarks/{type}/bam_rg_{sample}_{mapper}.benchmark.tsv"
    conda:
        "../Envs/bam_rg.yaml"
    threads: 4
    resources:
        mem_mb= 200000,
        runtime= 60,
        cpus_per_task= 4       
    shell:
        """
        picard AddOrReplaceReadGroups -XX:ActiveProcessorCount={threads} -Xmx10g \
            I={input.bamIn} \
            O={output.bamOut} \
            RGLB={params.sample} \
            RGPL=illumina \
            RGPU=unit \
            RGSM={params.sample} \
            RGID={params.sample} \
            2> {log}
        """


rule stats:
    params:
        mapper=MAPPER
    input:
        bamIn=expand("{output_dir}/Mapping/Bamout/{{mapper}}/{{type}}/mapping_rg_{{sample}}.bam",output_dir=config["output_dir"])
    output:
        statstsv=expand("{output_dir}/Analysis/{{mapper}}/{{type}}/perSample/{{sample}}.tsv",output_dir=config["output_dir"])
    log: 
        out= expand("{output_dir}/Logs/Analysis/{{type}}/stats_{{mapper}}_{{sample}}.out.log",output_dir=config["output_dir"]),
        err= expand("{output_dir}/Logs/Analysis/{{type}}/stats_{{mapper}}.{{sample}}.err.log",output_dir=config["output_dir"])
    benchmark:
       "../Benchmarks/{type}/stats_{mapper}_{sample}.benchmark.tsv"
    conda: 
        "../Envs/stats.yaml"
    #threads: NULL
    resources:
        mem_mb= 4000,
        runtime= 30,
        cpus_per_task= 1       
    shell:
        """
        python Scripts/Stats.py \
            -i {input.bamIn} \
            -o {output.statstsv} \
            > {log.out} \
            2> {log.err}
        """


rule merge_stats:
    params:
        mapper=MAPPER,
        inDirMonos=expand("{output_dir}/Analysis/{{mapper}}/monos/perSample/",output_dir=config["output_dir"]),
        inDirNonMonos=expand("{output_dir}/Analysis/{{mapper}}/nonmonos/perSample/",output_dir=config["output_dir"])

    input:
        statstsvMono=expand("{output_dir}/Analysis/{{mapper}}/{type}/perSample/{monos}.tsv",output_dir=config["output_dir"],monos=MONOS,type="monos"),
        statstsvNonMono=expand("{output_dir}/Analysis/{{mapper}}/{type}/perSample/{nonmonos}.tsv",output_dir=config["output_dir"],nonmonos=NONMONOS,type="nonmonos")

    output:
        statstsv=expand("{output_dir}/Analysis/{{mapper}}/stats.tsv",output_dir=config["output_dir"])
    conda: 
        "../Envs/statsCombine.yaml"
    #threads: NULL
    resources:
        mem_mb= 50000,
        runtime= 240,
        cpus_per_task= 1          
    shell:
        """
        Rscript Scripts/combineStatsFiles.R {params.inDirMonos} {output.statstsv} {params.inDirMonos} 
        """
rule filter:
    params: 
        mapper=MAPPER,
        filter_1=config["filter_1"],
        filter_2=config["filter_2"],
        filter_3=config["filter_3"],
        outprefix=expand("{output_dir}/Analysis/{{mapper}}/",output_dir=config["output_dir"])
    input: 
        statstsv=expand("{output_dir}/Analysis/{{mapper}}/stats.tsv",output_dir=config["output_dir"])
    output:
        Data1=expand("{output_dir}/Analysis/{{mapper}}/Data_1_Clusters_Target_vs_Reason_to_remove_8_15_1000_summed_per_species.txt", output_dir=config["output_dir"]),
        Data2=expand("{output_dir}/Analysis/{{mapper}}/Data_2_Clusters_filtered_due_to_homology_to_8_15_1000.txt", output_dir=config["output_dir"]),
        Data3=expand("{output_dir}/Analysis/{{mapper}}/Data_3_READ_COUNT_removed_CLUSTERS_SUM_8_15_1000.tsv", output_dir=config["output_dir"]),
        Data4=expand("{output_dir}/Analysis/{{mapper}}/Data_4_SUM_8_15_1000.tsv", output_dir=config["output_dir"]),
        Data5=expand("{output_dir}/Analysis/{{mapper}}/Data_5_SUM_MINREAD_FILTER_8_15_1000.tsv", output_dir=config["output_dir"])
    log: 
        out= expand("{output_dir}/Logs/Analysis/filter_{{mapper}}.out.log",output_dir=config["output_dir"]), 
        err= expand("{output_dir}/Logs/Analysis/filter_{{mapper}}.err.log",output_dir=config["output_dir"]),
        log= expand("{output_dir}/Logs/Analysis/filter_{{mapper}}.log.log",output_dir=config["output_dir"])
    benchmark:
       "../Benchmarks/filter_{mapper}.benchmark.tsv"
    conda: "../Envs/filter.yaml"
    resources:
        mem_mb= 50000,
        runtime= 10,
        cpus_per_task= 1       
    shell:
        """
        python Scripts/Parse_tsv.py \
            -i {input.statstsv} \
            -log {log.log} \
            -f1 {params.filter_1} \
            -f2 {params.filter_2} \
            -f3 {params.filter_3} \
            -op {params.outprefix} \
            > {log.out} \
            2> {log.err}
        """

# Finally, this rule attempts to summarize the results/execution success of a number of rules in this pipeline
# by using MultiQC over all generated logs. Granted, this skips over the majority of these logs, as 
# most of them are irrelevant, empty or perhaps only contain log-data originating from rules
# that generate logs that MultiQC is unable to parse.
# -----     
# Input:    - The files discussing or displaying specific data after applying user-defined filters to the tsv file.
#             [NOTE] These files are not actually used by this rule, and are only defined as input because 
#                    the existence of these files confirms the completion of the entire pipeline, which would in turn
#                    mean all log-files are available to parse at that moment.
# Output:   - A MultiQC output file in HTML format, summarizing all data it could parse from the entirity
#             of the log folder at the time of executing the command. 
rule logging:
    params:
        logdir=expand("{output_dir}/Logs/", output_dir=config["output_dir"])
    input:
        Data1=expand("{output_dir}/Analysis/{mapper}/Data_1_Clusters_Target_vs_Reason_to_remove_8_15_1000_summed_per_species.txt", output_dir=config["output_dir"], mapper=MAPPER),
        Data2=expand("{output_dir}/Analysis/{mapper}/Data_2_Clusters_filtered_due_to_homology_to_8_15_1000.txt", output_dir=config["output_dir"], mapper=MAPPER),
        Data3=expand("{output_dir}/Analysis/{mapper}/Data_3_READ_COUNT_removed_CLUSTERS_SUM_8_15_1000.tsv", output_dir=config["output_dir"], mapper=MAPPER),
        Data4=expand("{output_dir}/Analysis/{mapper}/Data_4_SUM_8_15_1000.tsv", output_dir=config["output_dir"], mapper=MAPPER),
        Data5=expand("{output_dir}/Analysis/{mapper}/Data_5_SUM_MINREAD_FILTER_8_15_1000.tsv", output_dir=config["output_dir"], mapper=MAPPER)
    output:
        alllogs=expand("{output_dir}/logSummary/multiQClogsummary.html", output_dir=config["output_dir"])
    log:
        out= expand("{output_dir}/Logs/Analysis/logging.out.log",output_dir=config["output_dir"]),
        err= expand("{output_dir}/Logs/Analysis/logging.err.log",output_dir=config["output_dir"])
    benchmark:
       "../Benchmarks/logging.benchmark.tsv"
    conda:
        "../Envs/logging.yaml"
    #threads: NULL
    resources:
        mem_mb= 1000,
        runtime= 10,
        cpus_per_task= 1          
    shell:
        """
        multiqc {params.logdir} \
            -n {output.alllogs} \
            > {log.out} \
            2> {log.err}
        """
