rule split_barcode_file:
    input:
        barcodes=expand("{input}/{barcodes}", input=config["input_dir"], barcodes=config["barcode_filename"])
    params:
        barcodesfiltered=expand("{output}/{barcodesfiltered}", output=config["output_dir"], barcodesfiltered=config["barcodefiltered_filename"])
    output:
        barcodesfiltered=expand("{output}/{barcodesfiltered}", output=config["output_dir"], barcodesfiltered=config["barcodefiltered_filename"])
    shell:    
        """
        rm -f {params.barcodesfiltered}
        header=$(awk 'NR==1 {{print; exit}}' {input.barcodes})
        content=$(awk 'NR=2 {{print; exit}}' {input.barcodes})
        IFS="	"; headerList=($header)

        for column in "${{!headerList[@]}}";
        do
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
        tail +2 {input.barcodes} | while read line; do
            if [ $headersPassed ]; then
                    R1=$( echo "$line" |cut -f $(($barcodeR1Index + 1)))
                    R2=$( echo "$line" |cut -f $(($barcodeR2Index + 1)))
                    Sample=$( echo "$line" |cut -f $(($sampleIndex + 1)))
                {{
                    printf "%s\t%s\t%s\n" "$R1" "$R2" "$Sample"
                }} >> {params.barcodesfiltered}
            else
                header=$( echo "$line")
                headersPassed=true
            fi
        done
        """   

rule deduplicate_trim:
    input:
        barcodes=expand("{input}/{barcodes}", input=config["input_dir"], barcodes=config["barcode_filename"]),
        barcodesfiltered=expand("{output}/{barcodesfiltered}",output=config["output_dir"], barcodesfiltered=config["barcodefiltered_filename"]),
        oriR1=expand("{input}/{read1}",input=config["input_dir"],read1=config["Read1"]),
        oriR2=expand("{input}/{read2}",input=config["input_dir"],read2=config["Read2"])
    params:
        read1=read1_sub,
        read2=read2_sub,
        barcodesfiltered=expand("{output}/{barcodesfiltered}",output=config["output_dir"], barcodesfiltered=config["barcodefiltered_filename"]),
        tmpdir=expand(tmpdirthis),
        outputdir=expand("{output}/Preprocessing", output=config["output_dir"]),
        inputdir=expand("{input}", input=config["input_dir"]),
        clonedir=expand("{tmp}/{clonefiltered}", tmp=tmpdirthis, clonefiltered=config["clonefiltered_dir"])
    output:
        filteredR1=temp(expand("{tmp}/{clonefiltered}/{name}.1.fq.gz", tmp=tmpdirthis, clonefiltered=config["clonefiltered_dir"], name=read1_sub)),
        filteredR2=temp(expand("{tmp}/{clonefiltered}/{name}.2.fq.gz", tmp=tmpdirthis, clonefiltered=config["clonefiltered_dir"], name=read2_sub)),
        demultiR1=temp(expand("{tmp}/{sample}.1.fq.gz", tmp=tmpdirthis, sample=SAMPLES)),
        demultiR2=temp(expand("{tmp}/{sample}.2.fq.gz", tmp=tmpdirthis, sample=SAMPLES))
    conda: "../Envs/deduplication.yaml"
    shell:
        """
        header=$(awk 'NR==1 {{print; exit}}' {input.barcodes})
        content=$(awk 'NR==2 {{print; exit}}' {input.barcodes})
        IFS="	"; headerList=($header)

        for column in "${{!headerList[@]}}";
        do
            case "${{headerList[$column]}}" in
                Wobble_R1)
                    WR1=$( echo "$content" |cut -f $(("${{column}}" + 1)))
                    ;;
                Wobble_R2)
                    WR2=$( echo "$content" |cut -f $(("${{column}}" + 1)))
                    ;;
                ENZ_R1)
                    ER1=$( echo "$content" |cut -f $(("${{column}}" + 1)))
                    ;;
                ENZ_R2)
                    ER2=$( echo "$content" |cut -f $(("${{column}}" + 1)))
                    ;;
            esac
        done
        fastp --in1 {input.oriR1} --in2 {input.oriR2} --out1 {params.clonedir}/{params.read1}.1.fq.gz --out2 {params.clonedir}/{params.read2}.2.fq.gz --adapter_sequence ATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT --adapter_sequence_r2 CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT --dedup --trim_poly_g --umi --umi_loc per_read --umi_len 3 
        process_radtags -1 {params.clonedir}/{params.read1}.1.fq.gz -2 {params.clonedir}/{params.read2}.2.fq.gz -b {params.barcodesfiltered} -o {params.tmpdir} -r -D --inline_inline --renz_1 "$ER1" --renz_2 "$ER2" --retain_header --disable_rad_check
        """

rule headers:
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
        headerR1=temp(expand("{tmp}/Preprocessing/Sampleheaders/{{sample}}.1.fq.gz",tmp=tmpdirthis)),
        headerR2=temp(expand("{tmp}/Preprocessing/Sampleheaders/{{sample}}.2.fq.gz",tmp=tmpdirthis))
    shell:
        """
        for i in 1 2; do
            zcat {params.Reads}."$i".fq.gz | cut -f1 -d ' ' | 
            sed -e '/^@/ s/$/\tRG:Z:{params.flowcell}_{params.lane}_{params.sample}/'| gzip -c   > {params.header}."$i".fq.gz
        done
        """

rule cat:
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
    shell:
        """
        for i in 1 2; do
            ls -v {params.path}Preprocessing/Sampleheaders/*."$i".fq.gz | xargs cat  >> {params.outbase}$i.fq.gz
        done
        """

#rule trim:
#    params:
#        sample='{sample}',
#        outputdir=expand("{path}",  path=tmpdirthis),
#        preprosup=expand("{preprosup}", preprosup=config["preprocessingsup_dir"])
#    input:
#        R1_in=expand("{path}/Preprocessing/Sampleheaders/{sample}.1.fq.gz",path=tmpdirthis,sample=SAMPLES),
#        R2_in=expand("{path}/Preprocessing/Sampleheaders/{sample}.2.fq.gz",path=tmpdirthis,sample=SAMPLES)
#    output:
#        R1_out=temp(expand("{path}/trimmed/{{sample}}.trimmed.pair1.truncated.gz",  path=tmpdirthis)),
#        R2_out=temp(expand("{path}/trimmed/{{sample}}.trimmed.pair2.truncated.gz",  path=tmpdirthis))
#    conda: "../Envs/read_trimming.yaml"
#    shell:
#        "AdapterRemoval --file1 {params.outputdir}/Preprocessing/Sampleheaders/{params.sample}.1.fq.gz "
#        "--file2 {params.outputdir}/Preprocessing/Sampleheaders/{params.sample}.2.fq.gz "
#        "--basename {params.outputdir}/trimmed/{params.sample}.trimmed --trimns --trimqualities "
#        "--adapter-list {params.preprosup}/adapters.txt --gzip"

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
    shell:
        """
        for i in 1 2; 
        do
            cp {params.sampleHeaders}."$i".fq.gz {params.outputdir}/{params.sample}.demultiplexed_R"$i".fq.gz
        done
        """
#mv {params.outputdir}."$i".fq.gz {params.outputdir}/{params.sample}.demultiplexed_R"$i".fq.gz

