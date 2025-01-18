# Authors v1.0 (Legacy):    ..., ..., ... & ...
# Authors v2.0 (New):       Jelle van der Heide
# Date:                     ../../..
# 
# This file contains rules for mapping with multiple different mappers. Currently, three options have 
# been defined. All three mappers require an indexation step prior to their execution. By adjusting 
# the [mapper_mode] option in the config-file, users can choose whether a specific mapper should be 
# used (Bowtie2, BWA or STAR) or if all three should be executed in the same run.
# --------------------------------------------------------------------------------------------------------------------


# This rule handles the indexing of the meta-reference file according to the way the Bowtie mapper expects
# it to be indexed. Additionally, as this is the start of the mapping process, a timestamp is written to a
# file for realtime benchmarking purposes. 
# -----     
# Input:    - The de-novo assembled meta-reference file.
# Output:   - An undisclosed number of index files, depending on reference size.
rule mapping_Bowtie2_index:
    params: 
        indexprefix=expand("{output_dir}/Mapping/Index/Bowtie/index", output_dir=config["output_dir"])
    input: 
        refblasted=ref
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
        32
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

# This rule executes a Bowtie mapping process, by mapping all preprocessed reads (this includes the
# reads that make up the meta-reference) on the de-novo assembled meta-reference. This requires 
# an undisclosed number of index files generated by the previous rule. 
# [NOTE] The current implementation is not set up for multimapping for benchmarking purposes.
#        To allow for multimapping, the parameter "-k 10" should be added to the shell block.
# -----
# Input:    - All index-files created by the previous rule.
#           - Both preprocessed read-files from one of the samples.
# Output:   - A sam-file containing an alignment for all reads from one of the samples on the meta-reference.
#           - A bam-file containing an alignment for all reads from one of the samples on the meta-reference.
#             [NOTE] this file does NOT have any readgroups just yet.
rule mapping_Bowtie2:   
    params:
        sample='{sample}',
        multimap=config["multimap_bowtie"],
        indexprefix=expand("{output_dir}/Mapping/Index/Bowtie/index", output_dir=config["output_dir"])
    input:
        index=expand("{output_dir}/Mapping/Index/Bowtie/index.1.bt2", output_dir=config["output_dir"]),
        r1=expand("{output_dir}/Preprocessing/Readgrouped/{{sample}}.1.fq.gz",output_dir=config["output_dir"]),
        r2=expand("{output_dir}/Preprocessing/Readgrouped/{{sample}}.1.fq.gz",output_dir=config["output_dir"])
    output:
        samOut=temp(expand("{tmp_dir}/Mapping/Samout/Bowtie/mapping_sq_{{sample}}.sam",tmp_dir=config["tmp_dir"])),
        bamOut=temp(expand("{tmp_dir}/Mapping/Bamout/Bowtie/mapping_sq_{{sample}}.bam",tmp_dir=config["tmp_dir"]))
    log: 
        bowtie2=expand("{output_dir}/Logs/Mapping/mapping_bowtie2_{{sample}}_bt.log",output_dir=config["output_dir"]),
        samtools=expand("{output_dir}/Logs/Mapping/mapping_bowtie2_{{sample}}_st.log",output_dir=config["output_dir"])
    benchmark:
       "../Benchmarks/mapping_Bowtie2_{sample}.benchmark.tsv"
    conda: 
        "../Envs/bowtie2.yaml"
    threads: 
        8
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


# This rule handles the indexing of the meta-reference file according to the way the Bwa mapper expects
# it to be indexed. Additionally, as this is the start of the mapping process, a timestamp is written to a
# file for realtime benchmarking purposes. 
# -----     
# Input:    - The de-novo assembled meta-reference file.
# Output:   - An undisclosed number of index files, depending on reference size.
rule mapping_bwa_index:
    params: 
        indexprefix=expand("{output_dir}/Mapping/Index/Bwa", output_dir=config["output_dir"])
    input: 
        refblasted=ref
    output:
        index=expand("{output_dir}/Mapping/Index/Bwa/index.amb", output_dir=config["output_dir"])
    log: 
        expand("{output_dir}/Logs/Mapping/mapping_bwa_index.log",output_dir=config["output_dir"])
    benchmark:
       "../Benchmarks/mapping_BWA_index.benchmark.tsv"
    conda: 
        "../Envs/bwa.yaml"
    threads: 
        32
    shell:
        """
        echo "Commencing Bwa mapping" >> time.txt
        date +%s%N >> time.txt
        bwa index \
            -p {params.indexprefix}/index \
            {input.refblasted} \
            2> {log}
        """

# This rule executes a Bwa mapping process, by mapping all preprocessed reads (this includes the
# reads that make up the meta-reference) on the de-novo assembled meta-reference. This requires 
# an undisclosed number of index files generated by the previous rule. 
# [NOTE] The current implementation is not set up for multimapping for benchmarking purposes.
#        To allow for multimapping, the parameter "-c 10" should be added to the shell block.
# -----         
# Input:    - All index-files created by the previous rule.
#           - Both preprocessed read-files from one of the samples.
# Output:   - A sam-file containing an alignment for all reads from one of the samples on the meta-reference.
#             A bam-file containing an alignment for all reads from one of the samples on the meta-reference.
#             [NOTE] this file does NOT have any readgroups just yet.
rule mapping_BWA:
    params:
        sample='{sample}',
        multimap=config["multimap_bwa"],
        indexprefix=expand("{output_dir}/Mapping/Index/Bwa/index", output_dir=config["output_dir"])
    input:
        index=expand("{output_dir}/Mapping/Index/Bwa/index.amb", output_dir=config["output_dir"]),
        r1=expand("{output_dir}/Preprocessing/Readgrouped/{{sample}}.1.fq.gz",output_dir=config["output_dir"]),
        r2=expand("{output_dir}/Preprocessing/Readgrouped/{{sample}}.1.fq.gz",output_dir=config["output_dir"])
    output:
        samOut=temp(expand("{tmp_dir}/Mapping/Samout/Bwa/mapping_sq_{{sample}}.sam",tmp_dir=config["tmp_dir"])),
        bamOut=temp(expand("{tmp_dir}/Mapping/Bamout/Bwa/mapping_sq_{{sample}}.bam",tmp_dir=config["tmp_dir"]))
    log: 
        bwa=expand("{output_dir}/Logs/Mapping/mapping_bwa_{{sample}}_bw.log",output_dir=config["output_dir"]),
        samtools=expand("{output_dir}/Logs/Mapping/mapping_bwa_{{sample}}_st.log",output_dir=config["output_dir"])
    benchmark:
       "../Benchmarks/mapping_BWA_{sample}.benchmark.tsv"
    conda: 
        "../Envs/bwa.yaml"
    threads: 
        8
    shell:
        """
        bwa mem \
            -t {threads} \
            -R '@RG\\tID:{params.sample}\\tSM:{params.sample}' \
            {params.indexprefix} \
            {params.multimap} \
            {input.r1} \
            {input.r2} \
            > {output.samOut} \
            2> {log.bwa}
        samtools view \
            -b \
            -o {output.bamOut} {output.samOut} \
            2> {log.samtools}
        """

# This rule handles the indexing of the meta-reference file according to the way the Star mapper expects
# it to be indexed. Additionally, as this is the start of the mapping process, a timestamp is written to a
# file for realtime benchmarking purposes. 
# -----     
# Input:    - The de-novo assembled meta-reference file.
# Output:   - An undisclosed number of index files, depending on reference size.
rule mapping_star_index:
    params:
        indexprefix=expand("{output_dir}/Mapping/Index/Star", output_dir=config["output_dir"]),
        indextemp=expand("../Misc/Mapping/Indexed/Star"),
        ram=config["star_ram"]
    input:
        refblasted=ref
    output:
        genome=expand("{output_dir}/Mapping/Index/Star/Genome" , output_dir=config["output_dir"])
    log: 
        expand("{output_dir}/Logs/Mapping/mapping_star_index.log",output_dir=config["output_dir"])
    benchmark: 
       "../Benchmarks/mapping_star_index.benchmark.tsv"
    conda: 
        "../Envs/star.yaml"
    threads: 
        32
    shell:
        """
        echo "Commencing STAR mapping" >> time.txt
        date +%s%N >> time.txt
        STAR \
            --genomeChrBinNbits 10 \
            --genomeSAindexNbases 10 \
            --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir {params.indexprefix} \
            --genomeFastaFiles {input.refblasted} \
            --outTmpDir {params.indextemp} \
            --limitGenomeGenerateRAM {params.ram} \
            2> {log}
        """

# This rule executes a Star mapping process, by mapping all preprocessed reads (this includes the
# reads that make up the meta-reference) on the de-novo assembled meta-reference. This requires 
# an undisclosed number of index files generated by the previous rule. 
# -----     
# Input:    - All index-files created by the previous rule.
#           - Both preprocessed read-files from one of the samples.
# Output:   - A sam-file containing an alignment for all reads from one of the samples on the meta-reference.
#             A bam-file containing an alignment for all reads from one of the samples on the meta-reference.
#             [NOTE] this file does NOT have any readgroups just yet.
rule mapping_Star:
    params:
        sample='{sample}',
        indexprefix=expand("{output_dir}/Mapping/Index/Star", output_dir=config["output_dir"]),
        r1out=expand("{output_dir}/Preprocessing/Readgrouped/{{sample}}.1.fq",output_dir=config["output_dir"]),
        r2out=expand("{output_dir}/Preprocessing/Readgrouped/{{sample}}.2.fq",output_dir=config["output_dir"]),
        logfinal=expand("{output_dir}/Mapping/{{sample)_Log.final", output_dir=config["output_dir"]),
        log=expand("{output_dir}/Mapping/{{sample)_Log", output_dir=config["output_dir"]),
        logprogress=expand("{output_dir}/Mapping/{{sample)_Log.progress", output_dir=config["output_dir"]),
        outtab=expand("{output_dir}/Mapping/{{sample)_SJ.out.tab", output_dir=config["output_dir"]),
        logout=expand("{output_dir}/Mapping/", output_dir=config["output_dir"])
    input:
        genome=expand("{output_dir}/Mapping/Index/Star/Genome" , output_dir=config["output_dir"]),
        r1=expand("{output_dir}/Preprocessing/Readgrouped/{{sample}}.1.fq.gz",output_dir=config["output_dir"]),
        r2=expand("{output_dir}/Preprocessing/Readgrouped/{{sample}}.2.fq.gz",output_dir=config["output_dir"])
    output:        
        samOut=temp(expand("{tmp_dir}/Mapping/Samout/Star/mapping_sq_{{sample}}.sam",tmp_dir=config["tmp_dir"])),
        bamOut=temp(expand("{tmp_dir}/Mapping/Bamout/Star/mapping_sq_{{sample}}.bam",tmp_dir=config["tmp_dir"]))
    log: 
        star=expand("{output_dir}/Logs/Mapping/mapping_star_{{sample}}_st.log",output_dir=config["output_dir"]),
        samtools=expand("{output_dir}/Logs/Mapping/mapping_star_{{sample}}_st.log",output_dir=config["output_dir"]),
    benchmark: 
       "../Benchmarks/mapping_star_{sample}.benchmark.tsv"
    conda: 
        "../Envs/star.yaml"
    threads: 
        8
    shell:
        """
        gunzip -f {input.r1}
        gunzip -f {input.r2}
        STAR \
            runThreadN {threads} \
            --genomeDir {params.indexprefix} \
            --readFilesIn {params.r1out} {params.r2out} \
            --outSAMattributes NM MD AS \
            --outSAMtype SAM \
            --outFileNamePrefix ../Output/Mapping/{params.sample}_ \
            --outFilterMatchNminOverLread 0.95 \
            --clip3pNbases 1 1 \
            --outSAMorder PairedKeepInputOrder \
            --outFilterMultimapScoreRange 0 \
            --alignEndsType Extend5pOfRead1 \
            --scoreGapNoncan 0 \
            --scoreGapGCAG 0 \
            --scoreGapATAC 0 \
            --scoreDelOpen 0 \
            --scoreDelBase 0 \
            --scoreInsOpen 0 \
            --scoreInsBase 0 \
            --alignMatesGapMax 20 \
            --readMapNumber \
            -1 2> {log.star}
        mv ../Output/Mapping/{params.sample}_Aligned.out.sam {output.samOut}
        mv {params.logfinal} {params.logout}
        mv {params.log} {params.logout}
        mv {params.logprogress} {params.logout}
        mv {params.outtab} {params.logout}
        samtools view -b -o {output.bamOut} {output.samOut} \
            2> {log.samtools}
        gzip -f {params.r1out}
        gzip -f {params.r2out}
        """
