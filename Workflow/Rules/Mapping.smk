rule mapping_Bowtie2_index:
    input: 
        refBlasted=expand("{path}/output_blast/Eukaryota_ref.fa",path=config["output_dir"]),
    output:
        temp("tmpfile.txt")
    conda: "../Envs/bowtie2.yaml"
    threads: 4
    shell:
        """
        bowtie2-build -f {input.refBlasted} ../Output/mapping/index -p 4
        touch tmpfile.txt
        """

rule mapping_Bowtie2:   
    params:
        sample='{sample}',
    input:
        mappingDone="tmpfile.txt",
        r1=expand("{tmp}/Preprocessing/Sampleheaders/{{sample}}.1.fq.gz",  tmp=tmpdirthis),#, sample=SAMPLES),
        r2=expand("{tmp}/Preprocessing/Sampleheaders/{{sample}}.2.fq.gz",  tmp=tmpdirthis)#, sample=SAMPLES)
    conda: "../Envs/bowtie2.yaml"
    output:
        samOut=temp(expand("{path}/mapping/mapping_sq_{{sample}}.sam",path=config["output_dir"])),
        bamOut=temp(expand("{path}/mapping/mapping_sq_{{sample}}.bam",path=config["output_dir"]))
    threads: 10
    shell:
        """
        bowtie2 -x ../Output/mapping/index -1 {input.r1} -2 {input.r2} -q --end-to-end --very-fast -k 10 --threads 10 -S {output.samOut}
        samtools view -b -o {output.bamOut} {output.samOut}
        """

rule mapping_BWA:

rule mapping_GEM3:

rule mapping_Segemehl:

