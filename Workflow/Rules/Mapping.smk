#rule mapping_Bowtie2_index:
#    input: 
#        refBlasted=expand("{path}/output_blast/Eukaryota_ref.fa",path=config["output_dir"]),
#    output:
#        index=expand("{path}/mapping/index.1.bt2", path=config["output_dir"])
#    log: "../Logs/Mapping/mapping_bowtie2_index.log"
#    conda: "../Envs/bowtie2.yaml"
#    threads: 4
#    shell:
#        """
#        echo "Commencing Bowtie mapping" >> time.txt
#        date +%s%N >> time.txt
#        bowtie2-build -f {input.refBlasted} ../Output/mapping/index -p 4 2> {log}
#        """ # Hier is nog enige output in de console die naar de log moet worden gestuurd

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
#        bamOut=expand("{path}/mapping/mapping_sq_{{sample}}.bam",path=config["output_dir"])
#    threads: 4
#    log: "../Logs/Mapping/mapping_bowtie2_{sample}.log"
#    shell:
#        """
#        bowtie2 -x ../Output/mapping/index -1 {input.r1} -2 {input.r2} -q --end-to-end --very-fast --threads 4 -S {output.samOut} 2> {log}
#        samtools view -b -o {output.bamOut} {output.samOut}
#        """ #k10

#rule mapping_bwa_index:
#    input: 
#        refBlasted=expand("{path}/output_blast/Eukaryota_ref.fa",path=config["output_dir"]),
#    output:
#        index=expand("{path}/mapping/index.amb", path=config["output_dir"])
#    conda: "../Envs/bwa.yaml"
#    log: "../Logs/Mapping/mapping_bwa_index.log"
#    threads: 4
#    shell:
#        """
#        echo "Commencing Bwa mapping" >> time.txt
#        date +%s%N >> time.txt
#        bwa index -p ../Output/mapping/index {input.refBlasted} 2> {log}
#        """

#rule mapping_BWA:
#    params:
#        sample='{sample}',
#    input:
#        index=expand("{path}/mapping/index.amb", path=config["output_dir"]),
#        r1=expand("{tmp}/Preprocessing/Sampleheaders/{{sample}}.1.fq.gz",  tmp=config["tmp_dir"]),#, sample=SAMPLES),
#        r2=expand("{tmp}/Preprocessing/Sampleheaders/{{sample}}.2.fq.gz",  tmp=config["tmp_dir"])#, sample=SAMPLES)
#    conda: "../Envs/bwa.yaml"
#    benchmark:"../Benchmarks/BWA-{sample}.benchmark.tsv"
#    log: "../Logs/Mapping/mapping_bwa_{sample}.log"
#    output:
#        samOut=temp(expand("{path}/mapping/mapping_rg_{{sample}}.sam",path=config["output_dir"])),
#        bamOut=temp(expand("{path}/mapping/mapping_rg_{{sample}}.bam",path=config["output_dir"])),
#    threads: 4
#    shell:
#        """
#        bwa mem -t 4 -c 10 -R '@RG\\tID:{params.sample}\\tSM:{params.sample}' ../Output/mapping/index {input.r1} {input.r2} > {output.samOut} 2> {log}
#        samtools view -b -o {output.bamOut} {output.samOut}
#        """ #c 10

#rule mapping_GEM3:

#rule mapping_Segemehl:

rule mapping_star_index:
    params:
        tempMaps=expand("../Misc/mapping/indexed")
    input:
        refBlasted=expand("{path}/output_blast/Eukaryota_ref.fa",path=config["output_dir"]),
    output:
        genome=expand("{output}/mapping/index/Genome" , output=config["output_dir"]),
    conda: "../Envs/star.yaml"
    log: "../Logs/Mapping/mapping_star_index.log"
    shell:
        """
        echo "Commencing Bowtie mapping" >> time.txt
        date +%s%N >> time.txt
        STAR --genomeSAindexNbases 10 --runThreadN 4 --runMode genomeGenerate --genomeDir ../Output/mapping/index --genomeFastaFiles {input.refBlasted} --outTmpDir {params.tempMaps} -limitGenomeGenerateRAM 1600000000000 2> {log}
        """

rule mapping_Star:
    params:
        r1out=temp(expand("{tmp}/Preprocessing/Sampleheaders/{{sample}}.1.fq", tmp=config["tmp_dir"])),
        r2out=temp(expand("{tmp}/Preprocessing/Sampleheaders/{{sample}}.2.fq", tmp=config["tmp_dir"])),
        sample='{sample}'
    input:
        genome=expand("{output}/mapping/index/Genome" , output=config["output_dir"]),
        refBlasted=expand("{path}/output_blast/Eukaryota_ref.fa",path=config["output_dir"]),
        r1=expand("{tmp}/Preprocessing/Sampleheaders/{{sample}}.1.fq.gz",  tmp=config["tmp_dir"]),#, sample=SAMPLES),
        r2=expand("{tmp}/Preprocessing/Sampleheaders/{{sample}}.2.fq.gz",  tmp=config["tmp_dir"])#, sample=SAMPLES)
    conda: "../Envs/star.yaml"
    log: "../Logs/Mapping/mapping_star_{sample}.log"
    output:        
        #r1out=temp(expand("{tmp}/Preprocessing/Sampleheaders/{{sample}}.1.fq", tmp=config["tmp_dir"])),
        #r2out=temp(expand("{tmp}/Preprocessing/Sampleheaders/{{sample}}.2.fq", tmp=config["tmp_dir"])),
        samOut=temp(expand("{path}/mapping/mapping_rg_{{sample}}.sam",path=config["output_dir"])),
        bamOut=expand("{path}/mapping/mapping_sq_{{sample}}.bam",path=config["output_dir"])
    threads: 16
    shell:
        """
	    gunzip {input.r1}
	    gunzip {input.r2}
        STAR runThreadN 16 --genomeDir ../Output/mapping/index --readFilesIn {params.r1out} {params.r2out} --outSAMattributes NM MD AS --outSAMtype SAM --outFileNamePrefix ../Output/mapping/{params.sample}_ --outFilterMatchNminOverLread 0.95 --clip3pNbases 1 1 --outSAMorder PairedKeepInputOrder --outFilterMultimapScoreRange 0 --alignEndsType Extend5pOfRead1 --scoreGapNoncan 0 --scoreGapGCAG 0 --scoreGapATAC 0 --scoreDelOpen 0 --scoreDelBase 0 --scoreInsOpen 0 --scoreInsBase 0 --alignMatesGapMax 20 --readMapNumber -1 2> {log}
        mv ../Output/mapping/{params.sample}_Aligned.out.sam {output.samOut}
        samtools view -b -o {output.bamOut} {output.samOut}
        """