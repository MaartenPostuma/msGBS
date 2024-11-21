# Authors v1.0 (Legacy): 
# Authors v2.0 (New):       Jelle van der Heide
# Date:                     ../../..
# 
# This file contains all msGBS pipeline functionality related to read pre-processing. 
# -------------------------------------------------------------------------------------------------------------------




# Considering the barcodefile contains a significant portion of data that is to be considered irrelevant,
# this rule handles extracting relevant data from said file. This data is stored in a new, split-up barcodefile which
# can in turn be used as input for other Rules instead of the bloated barcodefile.
# -----
# Input:    - A .tsv-file containing barcodes, which includes data surrounding their origin.
# Output:   - A filtered version of the file mentioned above. Data irrelevant to the pre-processing-steps
#             was removed from this version of the file to optimize data-extraction.
rule split_barcode_file:
    params:
        tmpdir=config["tmp_dir"]
    input:
        barcodefile=expand("{input}/{barcodes}", input=config["input_dir"], barcodes=config["barcode_filename"]),
        filteredR1=(expand("{tmp}/{clonefiltered}/{name}.1.fq.gz", tmp=config["tmp_dir"], clonefiltered=config["clonefiltered_dir"], name=read1_sub)),
        filteredR2=(expand("{tmp}/{clonefiltered}/{name}.2.fq.gz", tmp=config["tmp_dir"], clonefiltered=config["clonefiltered_dir"], name=read2_sub))
    output:
        barcodefilefiltered=temp(expand("{output}/{barcodesfiltered}", output=config["output_dir"], barcodesfiltered=config["barcodefiltered_filename"])),
        demultiR1=(expand("{tmp}/{sample}.1.fq.gz", tmp=config["tmp_dir"], sample=SAMPLES)),
        demultiR2=(expand("{tmp}/{sample}.2.fq.gz", tmp=config["tmp_dir"], sample=SAMPLES))
    log: "../Logs/Preprocessing/split_barcode_file.log"
    benchmark: "../Benchmarks/split_barcode_file.benchmark.tsv"
    conda: "../Envs/deduplication.yaml"
    threads: 8
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
        
        process_radtags -1 {input.filteredR1} -2 {input.filteredR2} -b {output.barcodefilefiltered} -o {params.tmpdir} -r -D --inline_inline --renz_1 "$ER1" --renz_2 "$ER2" --retain_header --disable_rad_check
        """   


# All reads have to be trimmed for adapters and Poly-G's. Besides trimming, PCR duplicates and low-quality reads
# are filtered out of the read files. with fastP. Afterwards, stacks is used to process radtags, which can be fetched 
# from the filtered barcode-file.
# -----
# Input:    - A .tsv-file containing barcodes, which includes data surrounding their origin.
#           - A filtered version of the file mentioned above. Data irrelevant to the pre-processing-steps.
#             was removed from this version of the file to optimize data-extraction.
#           - The first read file with raw reads.
#           - The second read file with raw reads.
# Output:   - The first file with trimmed reads, from which low quality reads and PCR duplicates have been removed. 
#           - The second file with trimmed reads, from which low quality reads and PCR duplicates have been removed. 
#           - Demultiplexed read files.
#           - 
rule deduplicate_trim:
    params:
        adapter1=config["adapter1"],
        adapter2=config["adapter2"]
    input:
        oriR1=expand("{input}/{read1}",input=config["input_dir"],read1=config["Read1"]),
        oriR2=expand("{input}/{read2}",input=config["input_dir"],read2=config["Read2"])
    output:
        filteredR1=(expand("{tmp}/{clonefiltered}/{name}.1.fq.gz", tmp=config["tmp_dir"], clonefiltered=config["clonefiltered_dir"], name=read1_sub)),
        filteredR2=(expand("{tmp}/{clonefiltered}/{name}.2.fq.gz", tmp=config["tmp_dir"], clonefiltered=config["clonefiltered_dir"], name=read2_sub))
    log: "../Logs/Preprocessing/deduplicate_trim.log"
    benchmark: "../Benchmarks/deduplicate_trim.benchmark.tsv"
    conda: "../Envs/deduplication.yaml"
    threads: 8
    shell:
        """
        fastp --in1 {input.oriR1} --in2 {input.oriR2} --out1 {output.filteredR1} --out2 {output.filteredR2} --adapter_sequence {params.adapter1} --adapter_sequence_r2 {params.adapter2} --dedup --trim_poly_g --umi --umi_loc per_read --umi_len 3 -j ../Output/Preprocessing.json -h ../Output/Preprocessing.html -w 8
        """ 


rule headers: 
    params:
        readfile="{readfile}",
        sample="{sample}",
        flowcell=flowCell,
        lane=lane,
        header=expand("{tmp}/Preprocessing/Sampleheaders/{{sample}}", tmp=config["tmp_dir"]),
        Reads=expand("{tmp}/{{sample}}", tmp=config["tmp_dir"], sample=SAMPLES),
        tmpdir=config["tmp_dir"]
    input:
        barcodesfiltered=expand("{output}/{barcodesfiltered}", output=config["output_dir"], barcodesfiltered=config["barcodefiltered_filename"]),
        Reads=expand("{tmp}/{sample}.{readfile}.fq.gz", tmp=config["tmp_dir"], sample=SAMPLES, readfile=readfile),
    output:
        Header=(expand("{tmp}/Preprocessing/Sampleheaders/{{sample}}.{{readfile}}.fq.gz",tmp=config["tmp_dir"])),
    log: "../Logs/Preprocessing/headers_{sample}_{readfile}.log"
    benchmark: "../Benchmarks/headers_{sample}_{readfile}.benchmark.tsv"
    #conda:
    threads: 1
    shell:
        """
        zcat {params.Reads}.{params.readfile}.fq.gz | cut -f1 -d ' ' | 
        sed -e '/^@/ s/$/\tRG:Z:{params.flowcell}_{params.lane}_{params.sample}/'| gzip -c   > {params.header}.{params.readfile}.fq.gz
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
        ls -v {params.path}Preprocessing/Sampleheaders/*.{readfile}.fq.gz | xargs cat  >> {params.outbase}{readfile}.fq.gz
        """

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
