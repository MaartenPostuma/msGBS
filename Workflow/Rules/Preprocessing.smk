# Authors v1.0 (Legacy): 
# Authors v2.0 (New):       Jelle van der Heide
# Date:                     ../../..
# 
# This file contains all msGBS pipeline functionality related to pre-processing raw reads in preparation for analysis. 
# --------------------------------------------------------------------------------------------------------------------



# All reads are trimmed for adapters and Poly-G's. Besides trimming, PCR duplicates and low-quality reads
# are filtered out of the read files with fastP.
# -----
# Input:    - The first read file with raw reads.
#           - The second read file with raw reads.
# Output:   - The first file with trimmed reads, from which low quality reads and PCR duplicates have been removed. 
#           - The second file with trimmed reads, from which low quality reads and PCR duplicates have been removed. 
#           - A report in HTML format containing quality control results.
#           - A report in JSON format containing quality control results.
rule deduplicate_trim:
    params:
        adapter1=config["adapter1"],
        adapter2=config["adapter2"]
    input:
        reads1=expand("{input}/{read1}",input=config["input_dir"],read1=config["Read1"]),
        reads2=expand("{input}/{read2}",input=config["input_dir"],read2=config["Read2"])
    output:
        filtered1=expand("{deduplicated_dir}/{name}.1.fq.gz", deduplicated_dir=config["deduplicated_dir"], name=read1_sub),
        filtered2=expand("{deduplicated_dir}/{name}.2.fq.gz", deduplicated_dir=config["deduplicated_dir"], name=read2_sub),
        fastp_html=expand("{preprocessing_out}/Preprocessing.html", preprocessing_out=config["preprocessing_out"]),
        fastp_json=expand("{preprocessing_out}/Preprocessing.json", preprocessing_out=config["preprocessing_out"]),
    log: 
        "../Logs/Preprocessing/deduplicate_trim.log"
    benchmark: 
        "../Benchmarks/deduplicate_trim.benchmark.tsv"
    conda: 
        "../Envs/deduplication.yaml"
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

# Considering the barcodefile contains a significant portion of data that is to be considered irrelevant,
# this rule handles extracting said relevant data from the file. This data is stored in a new, split-up barcodefile which
# can in turn be used as input for other Rules instead of the bloated original barcodefile. Afterwards, stacks is used to process radtags, which can be fetched 
# from the filtered barcode-file.
# -----     
# Input:    - A .tsv-file containing read barcodes, which includes data surrounding their origin.
#           - The first quality controlled read file.
#           - The second quality controlled read file.
# Output:   - A filtered version of the file mentioned above. Data irrelevant to the pre-processing-steps
#             was removed from this version of the file to optimize data-extraction.
#           - The first demultiplexed read file.
#           - The second demultiplexed read file.
rule split_barcodes_demultiplex:
    params:
        demultiplexed_dir=config["demultiplexed_dir"]
    input:
        barcodefile=expand("{input}/{barcodes}", input=config["input_dir"], barcodes=config["barcode_file"]),
        filtered1=expand("{deduplicated_dir}/{name}.1.fq.gz", deduplicated_dir=config["deduplicated_dir"], name=read1_sub),
        filtered2=expand("{deduplicated_dir}/{name}.2.fq.gz", deduplicated_dir=config["deduplicated_dir"], name=read2_sub),
    output:
        barcodefilefiltered=expand("{output}/{barcodesfiltered}", output=config["output_dir"], barcodesfiltered=config["barcodefiltered_file"]),
        demultiplexed1=expand("{demultiplexed_dir}/{sample}.1.fq.gz", demultiplexed_dir=config["demultiplexed_dir"], sample=SAMPLES),
        demultiplexed2=expand("{demultiplexed_dir}/{sample}.2.fq.gz", demultiplexed_dir=config["demultiplexed_dir"], sample=SAMPLES)
    log: 
        "../Logs/Preprocessing/split_barcode_file.log"
    benchmark: 
        "../Benchmarks/split_barcode_file.benchmark.tsv"
    conda: 
        "../Envs/deduplication.yaml"
    threads: 
        8
    shell:    
        """
        header=$(awk 'NR==1 {{print; exit}}' {input.barcodefile})
        content=$(awk 'NR==2 {{print; exit}}' {input.barcodefile})
        IFS="	"; headerList=($header)

        for column in "${{!headerList[@]}}"; do
            case "${{headerList[$column]}}" in
                Barcode_R1)
                    barcodeR1Index="${{column}}"
                    ;;
                Barcode_R2)
                    barcodeR2Index="${{column}}"
                    ;;
                Sample)
                    sampleIndex="${{column}}"
                    ;;
                ENZ_R1)
                    ER1=$( echo "$content" |cut -f $(("${{column}}" + 1)))
                    ;;
                ENZ_R2)
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
                {{
                    printf "%s\t%s\t%s\n" "$R1" "$R2" "$Sample"
                }} >> {output.barcodefilefiltered}
            else
                header=$( echo "$line")
                headersPassed=true
            fi
        done
        
        process_radtags \
            -1 {input.filtered1} \
            -2 {input.filtered2} \
            -b {output.barcodefilefiltered} \
            -o {params.demultiplexed_dir} \
            -r \
            --inline_inline \
            --renz_1 "$ER1" \
            --renz_2 "$ER2" \
            --retain_header \
            --disable_rad_check \
            --threads {threads} \
            2> {log}
        """

# To be able to properly track which segment of refseq the sample reads map on, the read headers have to be altered. This rule
# adjusts read headers so that it contains information like flowcell, lane and sample of origin as read group indicator.
# Considering there are multiple samples, as well as multiple versions of each sample read-file (as they are paired) 
# this rule repeats it's process for each version of each sample.
# -----     
# Input:    - A .tsv-file containing read barcodes, which includes data surrounding their origin.
#           - A quality controlled read file.
# Output:   - A filtered version of the file mentioned above. Data irrelevant to the pre-processing-steps
#             was removed from this version of the file to optimize data-extraction.
#           - The first demultiplexed read file.
rule readgroup_headers: 
    params:
        readfile="{readfile}",
        sample="{sample}",
        flowcell=flowCell,
        lane=lane,
    input:
        read=expand("{demultiplexed_dir}/{{sample}}.{{readfile}}.fq.gz", demultiplexed_dir=config["demultiplexed_dir"]),
    output:
        readgrouped=expand("{readgrouped_dir}/{{sample}}.{{readfile}}.fq.gz",readgrouped_dir=config["readgrouped_dir"]),
    #log: NULL
    benchmark: 
        "../Benchmarks/headers_{sample}_{readfile}.benchmark.tsv"
    #conda: NULL
    threads: 
        8
    shell:
        """
        zcat {input.read} | 
        cut -f1 -d ' ' | 
        sed -e '/^@/ s/$/\tRG:Z:{params.flowcell}_{params.lane}_{params.sample}/'| 
        gzip -c   > {output.readgrouped}
        """

rule cat:
    params:
        readfile="{readfile}",
        path=expand("{path}/", path=config["tmp_dir"]),
        outbase=expand("{path}/Preprocessing/demultiplexed_R",  path=config["output_dir"])
    input:
        R=expand("{path}/Preprocessing/Sampleheaders/{sample}.{readfile}.fq.gz",path=config["tmp_dir"], sample=SAMPLES, readfile=readfile),
    output:
        outFile=expand("{path}/Preprocessing/demultiplexed_R{{readfile}}.fq.gz",  path=config["output_dir"]),
    log: "../Logs/Preprocessing/cat_{readfile}.log"
    benchmark: "../Benchmarks/cat_{readfile}.benchmark.tsv"
    threads: 1
    shell:
        """
        ls -v {params.path}Preprocessing/Sampleheaders/*.{readfile}.fq.gz | 
        xargs cat  >> {params.outbase}{readfile}.fq.gz
        """

#rule cat:    backup
#    params:
#        readfile="{readfile}",
#        path=expand("{path}/", path=config["tmp_dir"]),
#        outbase=expand("{path}/Preprocessing/demultiplexed_R",  path=config["output_dir"])
#    input:
#        R=expand("{path}/Preprocessing/Sampleheaders/{sample}.{readfile}.fq.gz",path=config["tmp_dir"], sample=SAMPLES, readfile=readfile),
#    output:
#        outFile=expand("{path}/Preprocessing/demultiplexed_R{{readfile}}.fq.gz",  path=config["output_dir"]),
#    log: "../Logs/Preprocessing/cat_{readfile}.log"
#    benchmark: "../Benchmarks/cat_{readfile}.benchmark.tsv"
#    threads: 1
#    shell:
#        """
#        ls -v {params.path}Preprocessing/Sampleheaders/*.{readfile}.fq.gz | 
#        xargs cat  >> {params.outbase}{readfile}.fq.gz
#        """

rule move_monos: 
    params:
        readfile="{readfile}",
        sample='{sample}',
        outputdir=expand("{path}/Preprocessing/monos", path=config["output_dir"]),
    input:
        R=expand("{path}/Preprocessing/Sampleheaders/{{sample}}.{{readfile}}.fq.gz",path=config["tmp_dir"]),
    output:
        R_demulti=expand("{path}/Preprocessing/monos/{{sample}}.demultiplexed_R{{readfile}}.fq.gz",  path=config["output_dir"]),
    log: "../Logs/Preprocessing/move_monos_{sample}_{readfile}.log"
    benchmark: "../Benchmarks/move_monos_{sample}_{readfile}.benchmark.tsv"
    #conda:
    threads: 1
    shell:
        """
        cp {input.R} {output.R_demulti}
        """

