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
        cp {input.monoReads_R1} {output.monoRerads_R1}
        cp {input.monoReads_R2} {output.monoRerads_R2}
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
        cp {input.monoReads_R1} {output.monoRerads_R1}
        cp {input.monoReads_R2} {output.monoRerads_R2}
        """

rule createRef:
    input:
        indRef=expand("{ref_loc}/{monos}.fa",ref_loc=config["ref_loc"],monos=MONOS),
        checkMonos=expand("{output_dir}/Preprocessing/monoCheck.txt",output_dir=config["output_dir"]),
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
        runtime= 10,
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
        inDir=expand("{output_dir}/Analysis/{{mapper}}/perSample/",output_dir=config["output_dir"])
    input:
        statstsvMono=expand("{output_dir}/Analysis/{{mapper}}/{type}/perSample/{monos}.tsv",output_dir=config["output_dir"],sample=MONOS,type="monos"),
        statstsvNonMono=expand("{output_dir}/Analysis/{{mapper}}/{type}/perSample/{nonmonos}.tsv",output_dir=config["output_dir"],sample=NONMONOS,type="nonmonos")

    output:
        statstsv=expand("{output_dir}/Analysis/{{mapper}}/stats.tsv",output_dir=config["output_dir"])
    log: 
        out= expand("{output_dir}/Logs/Analysis/{{type}}/stats_{{mapper}}.out.log",output_dir=config["output_dir"]),
    benchmark:
       "../Benchmarks/{type}/stats_{mapper}.benchmark.tsv"
    conda: 
        "../Envs/statsCombine.yaml"
    #threads: NULL
    resources:
        mem_mb= 50000,
        runtime= 30,
        cpus_per_task= 1          
    shell:
        """
        echo 'Finished mapping with {params.mapper}'  >> time.txt 
        date +%s%N >> time.txt 
        Rscript Scripts/combineStatsFiles.R {params.inDir} {output.statstsv} > {log.out}
        """
