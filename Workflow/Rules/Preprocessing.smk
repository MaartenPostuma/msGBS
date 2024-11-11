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
    input:
        barcodefile=expand("{input}/{barcodes}", input=config["input_dir"], barcodes=config["barcode_filename"])
    output:
        barcodefilefiltered=temp(expand("{output}/{barcodesfiltered}", output=config["output_dir"], barcodesfiltered=config["barcodefiltered_filename"]))
    log: "../Logs/Preprocessing/split_barcode_file.log"
    benchmark: "../Benchmarks/split_barcode_file.benchmark.txt"
    #conda: None
    threads: 1
    shell:    
        """
        header=$(awk 'NR==1 {{print; exit}}' {input.barcodefile})
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
        tmpdir=expand(tmpdirthis),
    input:
        barcodes=expand("{input}/{barcodes}", input=config["input_dir"], barcodes=config["barcode_filename"]),
        barcodesfiltered=expand("{output}/{barcodesfiltered}",output=config["output_dir"], barcodesfiltered=config["barcodefiltered_filename"]),
        oriR1=expand("{input}/{read1}",input=config["input_dir"],read1=config["Read1"]),
        oriR2=expand("{input}/{read2}",input=config["input_dir"],read2=config["Read2"])
    output:
        filteredR1=(expand("{tmp}/{clonefiltered}/{name}.1.fq.gz", tmp=tmpdirthis, clonefiltered=config["clonefiltered_dir"], name=read1_sub)),
        filteredR2=(expand("{tmp}/{clonefiltered}/{name}.2.fq.gz", tmp=tmpdirthis, clonefiltered=config["clonefiltered_dir"], name=read2_sub)),
        demultiR1=(expand("{tmp}/{sample}.1.fq.gz", tmp=tmpdirthis, sample=SAMPLES)),
        demultiR2=(expand("{tmp}/{sample}.2.fq.gz", tmp=tmpdirthis, sample=SAMPLES))
    log: "../Logs/Preprocessing/deduplicate_trim.log"
    benchmark: "../Benchmarks/deduplicate_trim.benchmark.txt"
    conda: "../Envs/deduplication.yaml"
    threads: 1
    shell:
        """
        header=$(awk 'NR==1 {{print; exit}}' {input.barcodes})
        content=$(awk 'NR==2 {{print; exit}}' {input.barcodes})
        IFS="	"; headerList=($header)

        for column in "${{!headerList[@]}}";
        do
            case "${{headerList[$column]}}" in
                ENZ_R1)
                    ER1=$( echo "$content" |cut -f $(("${{column}}" + 1)))
                    ;;
                ENZ_R2)
                    ER2=$( echo "$content" |cut -f $(("${{column}}" + 1)))
                    ;;
            esac
        done
        fastp --in1 {input.oriR1} --in2 {input.oriR2} --out1 {output.filteredR1} --out2 {output.filteredR2} --adapter_sequence ATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT --adapter_sequence_r2 CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT --dedup --trim_poly_g --umi --umi_loc per_read --umi_len 3 -j ../Output/Preprocessing -h ../Output/Preprocessing
        process_radtags -1 {output.filteredR1} -2 {output.filteredR2} -b {input.barcodesfiltered} -o {params.tmpdir} -r -D --inline_inline --renz_1 "$ER1" --renz_2 "$ER2" --retain_header --disable_rad_check
        """ # radtags kan misschien naar de bovenstaande functie, dan kan ik die case samevoegen. is wat netter
            # op dit moment heb ik niet helemaal duidelijk of demultiplexing wel echt gebeurt hier.


rule headers: #parallelisatie
    params:
        sample="{sample}",
        flowcell=flowCell,
        lane=lane,
        header=expand("{tmp}/Preprocessing/Sampleheaders/{{sample}}", tmp=tmpdirthis),
        Reads=expand("{tmp}/{{sample}}", tmp=tmpdirthis, sample=SAMPLES),
        tmpdir=tmpdirthis
    input:
        barcodesfiltered=expand("{output}/{barcodesfiltered}", output=config["output_dir"], barcodesfiltered=config["barcodefiltered_filename"]),
        R1=expand("{tmp}/{sample}.1.fq.gz", tmp=tmpdirthis, sample=SAMPLES),
        R2=expand("{tmp}/{sample}.2.fq.gz", tmp=tmpdirthis, sample=SAMPLES)
    output:
        headerR1=(expand("{tmp}/Preprocessing/Sampleheaders/{{sample}}.1.fq.gz",tmp=tmpdirthis)),
        headerR2=(expand("{tmp}/Preprocessing/Sampleheaders/{{sample}}.2.fq.gz",tmp=tmpdirthis))
    log: "../Logs/Preprocessing/headers_{sample}.log"
    benchmark: "../Benchmarks/headers.benchmark.txt"
    #conda:
    threads: 1
    shell:
        """
        for i in 1 2; do
            zcat {params.Reads}."$i".fq.gz | cut -f1 -d ' ' | 
            sed -e '/^@/ s/$/\tRG:Z:{params.flowcell}_{params.lane}_{params.sample}/'| gzip -c   > {params.header}."$i".fq.gz
        done
        """

rule cat: # parallelisatie
    params:
        path=expand("{path}/", path=tmpdirthis),
        outbase=expand("{path}/Preprocessing/demultiplexed_R",  path=config["output_dir"])
    input:
        #R1_out=expand("{path}/trimmed/{sample}.trimmed.pair1.truncated.gz",path=tmpdirthis,sample=SAMPLES),
        #R2_out=expand("{path}/trimmed/{sample}.trimmed.pair2.truncated.gz",path=tmpdirthis,sample=SAMPLES)
        R1_in=expand("{path}/Preprocessing/Sampleheaders/{sample}.1.fq.gz",path=tmpdirthis, sample=SAMPLES),
        R2_in=expand("{path}/Preprocessing/Sampleheaders/{sample}.2.fq.gz",path=tmpdirthis, sample=SAMPLES)
    output:
        outFileR1=expand("{path}/Preprocessing/demultiplexed_R1.fq.gz",  path=config["output_dir"]),
        outFileR2=expand("{path}/Preprocessing/demultiplexed_R2.fq.gz",  path=config["output_dir"])
    log: "../Logs/Preprocessing/cat.log"
    benchmark: "../Benchmarks/cat.benchmark.txt"
    threads: 1
    shell:
        """
        for i in 1 2; do
            ls -v {params.path}Preprocessing/Sampleheaders/*."$i".fq.gz | xargs cat  >> {params.outbase}$i.fq.gz
        done
        """

rule move_monos: 
    params:
        sample='{sample}',
        inputdir=expand("{path}/trimmed",  path=tmpdirthis),
        outputdir=expand("{path}/Preprocessing/monos", path=config["output_dir"]),
        sampleHeaders=expand("{path}/Preprocessing/Sampleheaders/{{sample}}", path=tmpdirthis)
    input:
        #R1_in=expand("{path}/trimmed/{sample}.trimmed.pair1.truncated.gz",  path=tmpdirthis,sample=MONOS),
        #R2_in=expand("{path}/trimmed/{sample}.trimmed.pair2.truncated.gz",  path=tmpdirthis,sample=MONOS)
        R1_in=expand("{path}/Preprocessing/Sampleheaders/{{sample}}.1.fq.gz",path=tmpdirthis),
        R2_in=expand("{path}/Preprocessing/Sampleheaders/{{sample}}.2.fq.gz",path=tmpdirthis)
    output:
        R1_demulti=expand("{path}/Preprocessing/monos/{{sample}}.demultiplexed_R1.fq.gz",  path=config["output_dir"]),
        R2_demulti=expand("{path}/Preprocessing/monos/{{sample}}.demultiplexed_R2.fq.gz",  path=config["output_dir"])
    log: "../Logs/Preprocessing/move_monos_{sample}.log"
    benchmark: "../Benchmarks/move_monos.benchmark.txt"
    #conda:
    threads: 1
    shell:
        """
        for i in 1 2; 
        do
            cp {params.sampleHeaders}."$i".fq.gz {params.outputdir}/{params.sample}.demultiplexed_R"$i".fq.gz
        done
        """

