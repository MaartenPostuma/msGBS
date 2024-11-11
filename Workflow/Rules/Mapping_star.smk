#rule mapping_Star:
#    params:
#        tmpdirthis=tmpdirthis,
#        inputdir=expand("{path}/output_denovo",  path=config["output_dir"]),
#        outputdir=expand("{path}/mapping",  path=config["output_dir"]),
#        barcodes=expand("{path}/{bar}", path=config["input_dir"], bar=config["barcode_filename"]),
#    input:
#        refBlasted=expand("{path}/output_blast/Eukaryota_ref.fa",path=config["output_dir"]),
#        joined_All=expand("{path}/output_denovo/all.joined.fastq.gz",  path=config["output_dir"]),
#        merged_All=expand("{path}/output_denovo/all.merged.fastq.gz",  path=config["output_dir"])
#    output:
#        log=expand("{path}/mapping/mapping_variantcalling.log",path=config["output_dir"]),
#        bamOut=expand("{path}/mapping/out.bam",path=config["output_dir"])
#    threads: 40
#    shell:
#        "python Scripts/Map_STAR_snake.py --tmpdir {params.tmpdirthis} "
#        "--input_dir {params.inputdir} "
#        "--output_dir {params.outputdir} "
#        "--threads {threads} "
#        "--barcodes {params.barcodes}"


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


"""
rule merge:
    input:
        R1_in=expand("{path}/Preprocessing/demultiplexed_R1.fq.gz",  path=config["output_dir"]),
        R2_in=expand("{path}/Preprocessing/demultiplexed_R2.fq.gz",  path=config["output_dir"])
    output:
        unAssembled_1=temp(expand("{path}/output_denovo/unassembled_1.fastq.gz",  path=config["output_dir"])),
        unAssembled_2=temp(expand("{path}/output_denovo/unassembled_2.fastq.gz",  path=config["output_dir"])),
        merged_Assembled=expand("{path}/output_denovo/all.merged.fastq.gz",  path=config["output_dir"])
    params:
        unAssembled=expand("{path}/output_denovo/unassembled",  path=config["output_dir"]),
        merged=expand("{path}/output_denovo/all.merged.fastq.gz",  path=config["output_dir"]),
        preprosup=expand("{preprosup}", preprosup=config["preprocessingsup_dir"])

    threads: 32
    conda: "../Envs/merge_reads.yaml"
    shell:
        "NGmerge -1 {input.R1_in} -2 {input.R2_in} "
        "-o {params.merged} -n {threads} -f {params.unAssembled} "
        "-w {params.preprosup}/qual_profile.txt -q 33 -u 41 -z -g"
"""




"""
rule combine_joined_all:
    params:
        inputdir=expand("{path}/output_denovo/",  path=config["output_dir"]),
        cycles=getParam_cycle(param_cycle)
    input:
        unAssembled_1=expand("{path}/output_denovo/unassembled_1.fastq.gz",  path=config["output_dir"]),
        unAssembled_2=expand("{path}/output_denovo/unassembled_2.fastq.gz",  path=config["output_dir"]),
        barcodes=expand("{path}/{bar}", path=config["input_dir"], bar=config["barcode_filename"])
    output:
        joined_All=expand("{path}/output_denovo/all.joined.fastq.gz",  path=config["output_dir"])
    conda: "../Envs/combine_reads.yaml"
    shell:
        ""
        header=$(awk 'NR==1 {{print; exit}}' {input.barcodes})
        IFS="	"; headerList=($header)
        set +u
        for ((i=1; i<=${{#headerList[@]}}; i++)); do
            case "${{headerList[$i]}}" in
                Barcode_R1)
                    BR1i="$i"
                    ;;
                Barcode_R2)
                    BR2i="$i"
                    ;;
                Wobble_R1)
                    WR1i="$i"
                    ;;
                Wobble_R2)
                    WR2i="$i"
                    ;;
            esac
        done

        BR1MaxLen=0
        BR2MaxLen=0
        WR1Max=0
        WR2Max=0

        {{
            read
            while read line; do
                ifs="	"; thisLine=($line)
                if [ ${{thisLine[WR1i]}} -gt $WR1Max ]; then
                    WR1Max=${{thisLine[WR1i]}}
                fi
                if [ ${{thisLine[WR2i]}} -gt $WR2Max ]; then
                    WR2Max=${{thisLine[WR2i]}}
                fi
                if [ ${{#thisLine[BR1i]}} -gt $BR1MaxLen ]; then
                    BR1MaxLen=${{#thisLine[BR1i]}}
                fi
                if [ ${{#thisLine[BR2i]}} -gt $BR2MaxLen ]; then
                    BR2MaxLen=${{#thisLine[BR2i]}}
                fi
            done
        }} < {input.barcodes}
        maxR1=$((cycles - BR1MaxLen - WR1Max))
        maxR2=$((cycles - BR2MaxLen - WR2Max))
        paste <(seqtk seq {input.unAssembled_1} | cut -c1-$maxR1) <(seqtk seq -r {input.unAssembled_2} |cut -c1-$maxR2|seqtk seq -r -)|cut -f1-5|sed '/^@/!s/\t/NNNNNNNN/g'| sed s/+NNNNNNNN+/+/g| sed 's/ /\t/' | cut -f1,2 |  pigz -p 1 -c > {output.joined_All}
        ""
"""