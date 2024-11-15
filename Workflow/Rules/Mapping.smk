#rule mapping_Bowtie2_index:
#    input: 
#        refBlasted=expand("{path}/output_blast/Eukaryota_ref.fa",path=config["output_dir"]),
#    output:
#        index=expand("{path}/mapping/index.1.bt2", path=config["output_dir"])
#    conda: "../Envs/bowtie2.yaml"
#    threads: 4
#    shell:
#        """
#        bowtie2-build -f {input.refBlasted} ../Output/mapping/index -p 4
#        """

#rule mapping_Bowtie2:   
#    params:
#        sample='{sample}',
#    input:
#        index=expand("{path}/mapping/index.1.bt2",  path=config["output_dir"]),
#        r1=expand("{tmp}/Preprocessing/Sampleheaders/{{sample}}.1.fq.gz",  tmp=config["tmp_dir"]),#, sample=SAMPLES),
#        r2=expand("{tmp}/Preprocessing/Sampleheaders/{{sample}}.2.fq.gz",  tmp=config["tmp_dir"])#, sample=SAMPLES)
#    conda: "../Envs/bowtie2.yaml"
#    benchmark:"../Benchmarks/Bowtie2-{sample}.benchmark.tsv"
#    output:
#        samOut=temp(expand("{path}/mapping/mapping_sq_{{sample}}.sam",path=config["output_dir"])),
#        bamOut=temp(expand("{path}/mapping/mapping_sq_{{sample}}.bam",path=config["output_dir"]))
#    threads: 4
#    shell:
#        """
#        bowtie2 -x ../Output/mapping/index -1 {input.r1} -2 {input.r2} -q --end-to-end --very-fast --threads 4 -S {output.samOut}
#        samtools view -b -o {output.bamOut} {output.samOut}
#        """ #k10

rule mapping_bwa_index:
    input: 
        refBlasted=expand("{path}/output_blast/Eukaryota_ref.fa",path=config["output_dir"]),
    output:
        index=expand("{path}/mapping/index.amb", path=config["output_dir"])
    conda: "../Envs/bwa.yaml"
    threads: 4
    shell:
        """
        bwa index -p ../Output/mapping/index {input.refBlasted}
        touch tmpfile.txt
        """


rule mapping_BWA:
    params:
        sample='{sample}',
    input:
        index=expand("{path}/mapping/index.amb", path=config["output_dir"]),
        r1=expand("{tmp}/Preprocessing/Sampleheaders/{{sample}}.1.fq.gz",  tmp=config["tmp_dir"]),#, sample=SAMPLES),
        r2=expand("{tmp}/Preprocessing/Sampleheaders/{{sample}}.2.fq.gz",  tmp=config["tmp_dir"])#, sample=SAMPLES)
    conda: "../Envs/bwa.yaml"
    benchmark:"../Benchmarks/BWA-{sample}.benchmark.tsv"
    output:
        samOut=temp(expand("{path}/mapping/mapping_rg_{{sample}}.sam",path=config["output_dir"])),
        bamOut=temp(expand("{path}/mapping/mapping_rg_{{sample}}.bam",path=config["output_dir"])),
    threads: 4
    shell:
        """
        bwa mem -t 4 -c 10 -R '@RG\\tID:{params.sample}\\tSM:{params.sample}' ../Output/mapping/index {input.r1} {input.r2} > {output.samOut}
        samtools view -b -o {output.bamOut} {output.samOut}
        """ #c 10

#rule mapping_GEM3:

#rule mapping_Segemehl:
