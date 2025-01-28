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
    params:
        run="{run}",
        barcode_R1=config["barcode_R1"],
        barcode_R2=config["barcode_R2"],
        enzyme_R1=config["enzyme_R1"],
        enzyme_R2=config["enzyme_R2"],
        readfile_R1=config["readfile_R1"],
        readfile_R2=config["readfile_R2"],
        sample=config["sample"],
        readdir=directory(expand("{tmp_dir}/Preprocessing/Demultiplexed", tmp_dir=config["tmp_dir"]))
    input:
        barcodefile=expand("{input}/{barcodes}", input=config["input_dir"], barcodes=config["barcode_file"]),
    output:
        barcodefilefiltered=expand("{output_dir}/Preprocessing/Barcodesfiltered/{{run}}_barcodefiltered.tsv", output_dir=config["output_dir"]),
    log: 
        log1=expand("{output_dir}/Logs/Preprocessing/split_barcode_file_{{run}}.log", output_dir=config["output_dir"]),
        log2=expand("{output_dir}/Logs/Preprocessing/process_radtags_{{run}}.log", output_dir=config["output_dir"])
    benchmark: 
       "../Benchmarks/split_barcode_file.benchmark_{run}.tsv"
    resources:
        mem_mb= 1000,
        runtime= 2,
        cpus_per_task= 1
    conda:
        "../Envs/deduplication.yaml"
    threads: 
        1
    shell:    
        """
        header=$(awk 'NR==1 {{print; exit}}' {input.barcodefile})
        content=$(awk 'NR==2 {{print; exit}}' {input.barcodefile})
        IFS="	"; headerList=($header)

        for column in "${{!headerList[@]}}"; do
            case "${{headerList[$column]}}" in
                {params.barcode_R1})
                    barcodeR1Index="${{column}}"
                    ;;
                {params.barcode_R2})
                    barcodeR2Index="${{column}}"
                    ;;
                {params.sample})
                    sampleIndex="${{column}}"
                    ;;
                {params.readfile_R1})
                    rawR1Index="${{column}}"
                    ;;
                {params.enzyme_R1})
                    ER1=$( echo "$content" |cut -f $(("${{column}}" + 1)))
                    ;;
                {params.enzyme_R2})
                    ER2=$( echo "$content" |cut -f $(("${{column}}" + 1)))
                    ;;
            esac
        done
        headersPassed=false
        tail +2 {input.barcodefile} | while read line; do
            if [ $headersPassed ]; then
                    R1=$( echo "$line" |cut -f $(($barcodeR1Index + 1)))
                    R2=$( echo "$line" |cut -f $(($barcodeR2Index + 1)))
                    Sample=$( echo "$line" |cut -f $(($sampleIndex + 1)))
                    rawR1=$( echo "$line" |cut -f $(($rawR1Index + 1)))
                {{
                    if [[ $rawR1 == "{params.run}_R1.fq.gz" ]]; then
                        printf "%s\t%s\t%s\n" "$R1" "$R2" "$Sample"
                    fi
                }} >> {output.barcodefilefiltered}
            else
                header=$( echo "$line")
                headersPassed=true
            fi
        done
        """



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
        runtime= 120,
        cpus_per_task= 8
    threads: 
        8
    shell:    
        """
        mkdir {output.tmpdir}
        process_radtags \
            -1 {input.filtered1} \
            -2 {input.filtered2} \
            -b {output.barcodefilefiltered} \
            -o {output.tmpdir} \
            -r \
            --inline_inline \
            --renz_1 PacI \
            --renz_2 NsiI \
            --retain_header \
            --disable_rad_check \
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
        lambda w: {SAMPLES[w.sample]}, 
        readfile="{readfile}",
        tmp_dir=config["tmp_dir"],
        renamedunzipped=expand("{tmp_dir}/Preprocessing/Renamed/{{sample}}.{{readfile}}.fq", tmp_dir=config["tmp_dir"], sample=SAMPLES, readfile=readfile)
    input:
        tmpdir=expand("{tmp_dir}/Preprocessing/Demultiplexed/{run}", tmp_dir=config["tmp_dir"], run=RUN)
    output:
        renamed=temp(expand("{tmp_dir}/Preprocessing/Renamed/{{sample}}.{{readfile}}.fq.gz", tmp_dir=config["tmp_dir"], sample=SAMPLES, readfile=readfile))
    log:
        expand("{output_dir}/Logs/Preprocessing/rename_samples_{{sample}}_{{readfile}}.log", output_dir=config["output_dir"])
    benchmark: 
       "../Benchmarks/rename_samples_{sample}_{readfile}.benchmark.tsv"
    #conda: NULL
    resources:
        mem_mb= 1000,
        runtime= 5,
        cpus_per_task= 1
    threads:
        1
    shell:
        """
        zcat {params.tmp_dir}/Preprocessing/Demultiplexed/*/{wildcards.sample}.{params.readfile}.fq.gz > {params.renamedunzipped} 2> {log}
        gzip {params.renamedunzipped} 2> {log}
        """

# To be able to properly track which segment of refseq the sample reads map on, the read headers have to be altered. This rule
# adjusts read headers so that it contains information like flowcell, lane and sample of origin as read group indicator.
# -----     
# Input:    - A demultiplexed read file.
# Output:   - A demultiplexed read file readgroups added to it's read headers.
rule readgroup_headers: 
    params:
        sample="{sample}",
        readfile="{readfile}",
        flowcell=flowCell,
        lane=lane
    input:
        read=expand("{tmp_dir}/Preprocessing/Renamed/{{sample}}.{{readfile}}.fq.gz",tmp_dir=config["tmp_dir"])
    output:
        readgrouped=expand("{output_dir}/Preprocessing/Readgrouped/{{sample}}.{{readfile}}.fq.gz",output_dir=config["output_dir"])
    #log: NULL
    benchmark: 
       "../Benchmarks/headers_{sample}_{readfile}.benchmark.tsv"
    #conda: NULL
    resources:
        mem_mb= 1000,
        runtime= 10,
        cpus_per_task= 1
    threads: 
        1
    shell:
        """
        if [ -f '{input.read}' ]; then
            zcat {input.read} | 
            cut -f1 -d ' ' | 
            sed -e '/^@/ s/$/\tRG:Z:{params.flowcell}_{params.lane}_{params.sample}/'| 
            gzip -c   > {output.readgrouped}
        else
            touch {output.readgrouped}
        fi
        """

# This rule re-combines all reads into two fully preprocessed read files. While these files are currently unused 
# downstream, they are still generated to give users insight in how their input has been processed. 
# -----     
# Input:    - All preprocessed read files (demultiplexed), all either numbered 1 or 2.
# Output:   - A file which contains all preprocessed reads that remain, numbered either 1 or 2.
rule cat_all_preprocessed:
    params:
        readfile="{readfile}",
        tmp_dir=expand("{tmp_dir}/Preprocessing/Readgrouped/", tmp_dir=config["tmp_dir"])
    input:
        readgroupedsample=expand("{tmp_dir}/Preprocessing/Readgrouped/{sample}.{{readfile}}.fq.gz",tmp_dir=config["tmp_dir"], sample=SAMPLES)
    output:
        allpreprocessed=expand("{output_dir}/Output/Preprocessing/Preprocessed/preprocessed_R{{readfile}}.fq.gz",  output_dir=config["output_dir"])
    log: 
        expand("{output_dir}/Logs/Preprocessing/cat_{{readfile}}.log", output_dir=config["output_dir"])
    benchmark: 
       "../Benchmarks/cat_{readfile}.benchmark.tsv"
    #conda: NULL
    resources:
        mem_mb= 1000,
        runtime= 10,
        cpus_per_task= 1
    threads: 
        1
    shell:
        """
        ls -v {params.tmp_dir}*.{params.readfile}.fq.gz | 
        xargs cat  >> {output.allpreprocessed}
        """

# Since the reference-sequence is to be constructed exclusively from reads that come from mono-samples, these sample-reads are copied to another
# directory, where they can easily be called by reference creation rules.
# -----     
# Input:    - A preprocessed mono read-file.
# Output:   - A preprocessed mono read-file, now copied to a different directory.
rule move_monos: 
    input:
        mono=expand("{output_dir}/Preprocessing/Readgrouped/{{sample}}.{{readfile}}.fq.gz",output_dir=config["output_dir"])
    output:
        movedmono=expand("{output_dir}/Preprocessing/Preprocessedmonos/{{sample}}.{{readfile}}.fq.gz", output_dir=config["output_dir"])
    log: 
        expand("{output_dir}/Logs/Preprocessing/move_monos_{{sample}}_{{readfile}}.log", output_dir=config["output_dir"])
    benchmark: 
       "../Benchmarks/move_monos_{sample}_{readfile}.benchmark.tsv"
    #conda: NULL
    resources:
        mem_mb= 1000,
        runtime= 10,
        cpus_per_task= 1
    threads: 
        1
    shell:
        """
        cp {input.mono} {output.movedmono}
        """
