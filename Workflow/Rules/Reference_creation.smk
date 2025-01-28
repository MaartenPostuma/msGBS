# Authors v1.0 (Legacy):    ..., ..., ... & ...
# Authors v2.0 (New):       Jelle van der Heide
# Date:                     ../../..
# 
# To be able to determine the species of origin for the sample reads, reference sequences are required. However, as the majority of 
# refseqs for the target species are inadequate and much longer than needed the pipeline agenerates all de-novo from
# the mono samples. As for each mono there are two read files, they need to be merged or joined, depending on whether there is
# overlap between the reads. All mono reference is then merged into one large meta-reference file. 
# --------------------------------------------------------------------------------------------------------------------



# This rule merges paired mono-sample reads when there is enough overlap between them to do so. Reads that cannot be merged
# are placed in a different file within which they will end up being joined together instead.
# -----     
# Input:    - The preprocessed read files for one of the monos.
# Output:   - The unmerged reads from the preprocessed read files from one of the monos.
#           - The merged reads from the preprocessed read files from one of the monos.
rule merge_monos:
    params:
        supplement="Supplement",
        outputprefix=expand("{tmp_dir}/Reference_creation/Joined/{{sample}}.notmerged",  tmp_dir=config["tmp_dir"])
    input:
        R1=expand("{output_dir}/Preprocessing/Preprocessedmonos/{{sample}}.1.fq.gz",  output_dir=config["output_dir"]),
        R2=expand("{output_dir}/Preprocessing/Preprocessedmonos/{{sample}}.2.fq.gz",  output_dir=config["output_dir"])
    output:
        notmergedR1=temp(expand("{tmp_dir}/Reference_creation/Joined/{{sample}}.notmerged_1.fastq.gz",  tmp_dir=config["tmp_dir"])),
        notmergedR2=temp(expand("{tmp_dir}/Reference_creation/Joined/{{sample}}.notmerged_2.fastq.gz",  tmp_dir=config["tmp_dir"])),
        merged=temp(expand("{tmp_dir}/Reference_creation/Merged/{{sample}}.merged.fastq.gz",  tmp_dir=config["tmp_dir"]))
    log: 
        expand("{output_dir}/Logs/Reference_creation/merge_monos_{{sample}}.log",  output_dir=config["output_dir"])
    benchmark:
       "../Benchmarks/merge_monos_{sample}.benchmark.tsv"
    conda: 
        "../Envs/merge_reads.yaml"
    threads: 
        4
    resources:
        mem_mb= 10000,
        runtime= 60,
        cpus_per_task= 4
    shell:
        """
        NGmerge \
            -1 {input.R1} \
            -2 {input.R2} \
            -o {output.merged} \
            -n {threads} \
            -f {params.outputprefix} \
            -w {params.supplement}/qual_profile.txt \
            -q 33 \
            -u 41 \
            -z \
            -g \
            2> {log}
        """
        

# This rule joins together paired mono-sample reads when they could not be merged because of a lack of overlap between them. When 
# reads are joined together, eight 'N' characters are placed in between them, so them being in separate to an unknown extent is clear
# when mapping on these reads later on. To properly be able to join, some calculations must be made over information found in the barcode-file,
# as well as the number of cycles that have occurred during PCR.
# -----     
# Input:    - The unmerged reads from a mono-sample readfile.
#           - The barcode file for this run.
# Output:   - The joined reads from a mono-sample readfile.
rule join_monos:
    params:
        barcode_R1=config["barcode_R1"],
        barcode_R2=config["barcode_R2"],
        wobble_R1=config["wobble_R1"],
        wobble_R2=config["wobble_R2"],
        cycles=config["pcr_cycles"]
    input:
        notmergedR1=expand("{tmp_dir}/Reference_creation/Joined/{{sample}}.notmerged_1.fastq.gz",  tmp_dir=config["tmp_dir"]),
        notmergedR2=expand("{tmp_dir}/Reference_creation/Joined/{{sample}}.notmerged_2.fastq.gz",  tmp_dir=config["tmp_dir"]),
        barcodes=expand("{input_dir}/{barcodes}", input_dir=config["input_dir"], barcodes=config["barcode_file"])
    output:
        joined=temp(expand("{tmp_dir}/Reference_creation/Joined/{{sample}}.joined.fastq.gz",  tmp_dir=config["tmp_dir"]))
    #log: NULL
    benchmark:
       "../Benchmarks/join_monos_{sample}.benchmark.tsv"
    conda: 
        "../Envs/combine_reads.yaml"
    resources:
        mem_mb= 1000,
        runtime= 10,
        cpus_per_task= 1
    threads:
        1
    shell:
        """
        header=$(awk 'NR==1 {{print; exit}}' {input.barcodes})
        IFS="	"; headerList=($header)
        set +u
        for ((i=1; i<=${{#headerList[@]}}; i++)); do
            case "${{headerList[$i]}}" in
                {params.barcode_R1})
                    BR1i="$i"
                    ;;
                {params.barcode_R2})
                    BR2i="$i"
                    ;;
                {params.wobble_R1})
                    WR1i="$i"
                    ;;
                {params.wobble_R2})
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

        maxR1=$(({params.cycles} - BR1MaxLen - WR1Max))
        maxR2=$(({params.cycles} - BR2MaxLen - WR2Max))
        paste <(seqtk seq {input.notmergedR1} | cut -c1-141) <(seqtk seq -r {input.notmergedR2} | cut -c1-141 | seqtk seq -r -) | 
        cut -f1-5 | 
        sed $'/^@/!s/\t/NNNNNNNN/g' | 
        sed s/+NNNNNNNN+/+/g | 
        sed $'s/ /\t/' | 
        cut -f1,2 | 
        pigz -p 1 -c > {output.joined}
        """

# To be able to process the mono read-files, this rule combines the merged and joined reads into a single file for each mono.
# -----     
# Input:    - The joined reads for one of the monos.
#           - The merged reads for one of the monos.
# Output:   - A read file that contains all joined and merged reads for one of the monos.
rule cat_monos:
    input:
        merged=expand("{tmp_dir}/Reference_creation/Merged/{{sample}}.merged.fastq.gz",  tmp_dir=config["tmp_dir"]),
        joined=expand("{tmp_dir}/Reference_creation/Joined/{{sample}}.joined.fastq.gz",  tmp_dir=config["tmp_dir"])
    output:
        combined=temp(expand("{tmp_dir}/Reference_creation/Joined/{{sample}}.combined.fastq.gz",  tmp_dir=config["tmp_dir"]))
    #log: NULL
    benchmark:
       "../Benchmarks/join_monos_{sample}.benchmark.tsv"
    #conda: NULL
    resources:
        mem_mb= 1000,
        runtime= 10,
        cpus_per_task= 1
    #threads: NULL
    shell:
        """
        cat {input.joined} {input.merged} > {output.combined}
        """

# To prepare the reads for dereplication, they have to be sorted by length.
# -----     
# Input:    - The combined mono-reads.
# Output:   - A sorted version of the combined mono-reads.
rule sort_monos:
    input:
        combined=expand("{tmp_dir}/Reference_creation/Joined/{{sample}}.combined.fastq.gz",  tmp_dir=config["tmp_dir"])
    output:
        sorted=temp(expand("{tmp_dir}/Reference_creation/Sorted/{{sample}}.sorted.fa",  tmp_dir=config["tmp_dir"]))
    log: 
        expand("{output_dir}/Logs/Reference_creation/sort_monos_{{sample}}.log",  output_dir=config["output_dir"])
    benchmark:
       "../Benchmarks/sort_monos_{sample}.benchmark.tsv"
    conda: 
        "../Envs/sort.yaml"
    resources:
        mem_mb= 10000,
        runtime= 20,
        cpus_per_task= 4
    threads:
        4
    shell:
        """
        vsearch \
            -sortbylength {input.combined} \
            --output {output.sorted} \
            --threads {threads} \
            2> {log}
        """

# To clear out any duplicate reads before assembling the reference sequence, this rule dereplicates the
# reads. To be able to do this properly, it utilizes the sorted reads from the previous rule.
# -----     
# Input:    - The sorted mono-read files.
# Output:   - A dereplicated version of the sorted mono-read-files.
rule derep_monos:
    params:
        min_unique_size=config["min_unique_size"]
    input:
        sorted=expand("{tmp_dir}/Reference_creation/Sorted/{{sample}}.sorted.fa",  tmp_dir=config["tmp_dir"])
    output:
        dereplicated=temp(expand("{tmp_dir}/Reference_creation/Dereplicated/{{sample}}.dereplicated.fa",  tmp_dir=config["tmp_dir"]))
    log: 
        expand("{output_dir}/Logs/Reference_creation/derep_monos_{{sample}}.log",  output_dir=config["output_dir"])
    benchmark:
       "../Benchmarks/derep_monos_{sample}.benchmark.tsv"
    conda: 
        "../Envs/derep_monos.yaml"
    resources:
        mem_mb= 10000,
        runtime= 20,
        cpus_per_task= 4
    threads:
        4
    shell:
        """
        vsearch \
            -derep_fulllength {input.sorted} \
            -sizeout \
            -minuniquesize {params.min_unique_size} \
            -output {output.dereplicated} \
            --threads {threads} \
            2> {log}
        """

# To prepare the reads for clustering, they once again have to be sorted by length.
# -----     
# Input:    - A dereplicated version of the sorted mono-read-files.
# Output:   - A sorted version of the dereplicated mono-read files.
rule resort_monos:
    params:
        sample='{sample}'
    input:
        dereplicated=expand("{tmp_dir}/Reference_creation/Dereplicated/{{sample}}.dereplicated.fa", tmp_dir=config["tmp_dir"])
    output:
        sorted=temp(expand("{tmp_dir}/Reference_creation/Sorted/{{sample}}.resorted.fa", tmp_dir=config["tmp_dir"]))
    log: 
        expand("{output_dir}/Logs/Reference_creation/resort_monos_{{sample}}.log",  output_dir=config["output_dir"])
    benchmark:
       "../Benchmarks/resort_monos_{sample}.benchmark.tsv"
    conda: 
        "../Envs/sort.yaml"
    resources:
        mem_mb= 10000,
        runtime= 20,
        cpus_per_task= 4
    threads:
        4
    shell:
        """
        vsearch \
            -sortbylength {input.dereplicated} \
            --output {output.sorted} \
            --threads {threads} \
            2> {log}
        """

# In an effort to minimize duplicate reads before creating the de-novo reference, the reads are clustered.
# this allows the pipeline to omit reads that are too similar, as it is likely they are duplicate and might
# interfere with the mapping process.
# -----     
# Input:    - A sorted version of the dereplicated mono-read files.
# Output:   - A clustered version of the dereplicated mono-read files.
rule cluster:
    input:
        sorted=expand("{tmp_dir}/Reference_creation/Sorted/{{sample}}.resorted.fa", tmp_dir=config["tmp_dir"])
    output:
        clustered=temp(expand("{tmp_dir}/Reference_creation/Clustered/{{sample}}.clustered.fa", tmp_dir=config["tmp_dir"]))
    log: 
        expand("{output_dir}/Logs/Reference_creation/cluster_{{sample}}.log",  output_dir=config["output_dir"])
    benchmark:
       "../Benchmarks/cluster_{sample}.benchmark.tsv"
    conda: 
        "../Envs/cluster.yaml"
    resources:
        mem_mb= 10000,
        runtime= 20,
        cpus_per_task= 4
    threads: 
        4
    shell:
         """
         vsearch \
            -cluster_smallmem {input.sorted} \
            -id 0.95 \
            -centroids {output.clustered} \
            -sizeout \
            -strand both \
            --threads {threads} \
            2> {log}
        """

# This rule  adjusts read headers in preparation for assembling the reference sequence in a way that
# they only contain the origin and an index. Without this edit, the mapping results would be convoluted and 
# more complex to parse when analyzing it's results.
# -----     
# Input:    - A clustered version of the dereplicated mono-read files.
# Output:   - A clustered version of the dereplicated mono-read files, however here all headers are adjusted to 
#             contain it's origin and an index.
rule rename_fast:
    params:
        sample='{sample}'
    input:
        clustered=expand("{tmp_dir}/Reference_creation/Clustered/{{sample}}.clustered.fa", tmp_dir=config["tmp_dir"])
    output:
        renamed=temp(expand("{tmp_dir}/Reference_creation/Renamed/{{sample}}.renamed.fa", tmp_dir=config["tmp_dir"]))
    #log: NULL
    benchmark:
       "../Benchmarks/rename_fast_{sample}.benchmark.tsv"
    #conda: NULL
    resources:
        mem_mb= 10000,
        runtime= 20,
        cpus_per_task= 1
    #threads: NULL
    shell: 
        """
        index=0
        if [ {params.sample} = "chr" ]; then
            genotype={params.sample}
        elif [ {params.sample} ]; then
            genotype="{params.sample}""_"
        else 
            genotype=""
        fi
        while read line; do
            if [[ $line == \>* ]]; then
                index=$((index+1))
                echo -e ">$genotype$index" >> {output.renamed}
            elif [[ "$line" =~ ^@.* ]]; then
                if [ !$index ]; then
                    index=$((index+1))
                fi

                if [[ "$line" =~ "^@.*/[12]?\n$" ]]; then
                    end=${{line:-2}} 
                else 
                    end=""
                fi
                echo -e "@$genotype$index$end" >> {output.renamed}
            elif [[ "$line" =~ ^+.* ]]; then
                echo -e "+$genotype$index$end" >> {output.renamed}
            else
                echo "$line" >> {output.renamed}
            fi
        done < {input.clustered}
        """

# With this rule the de-novo meta-refseq is essentially created, as it combines all remaining dereplicated mono-reads
# into a single meta-reference file. 
# -----     
# Input:    - A clustered version of the dereplicated mono-read files, however here all headers are adjusted to 
#             contain it's origin and an index.
# Output:   - A de-novo meta reference file containing all dereplicated & clustered mono-reads.
rule ref_out:
    params:
        inputprefix=expand("{tmp_dir}/Reference_creation/Renamed/",  tmp_dir=config["tmp_dir"])
    input:
        renamed=expand("{tmp_dir}/Reference_creation/Renamed/{sample}.renamed.fa", tmp_dir=config["tmp_dir"], sample=MONOS)
    output:
        ref=expand("{output_dir}/Reference_creation/ref.fa", output_dir=config["output_dir"])
    #log: NULL
    benchmark:
       "../Benchmarks/ref_out.benchmark.tsv"
    resources:
        mem_mb= 10000,
        runtime= 20,
        cpus_per_task= 1
    shell:
        """
        cat {params.inputprefix}*.renamed.fa > {output.ref}
        """
