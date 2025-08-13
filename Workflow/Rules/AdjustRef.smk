rule move_monos_preprocess:
    input:
        monoReads_R1=expand("{reads_loc}/{{monos}}.1.fq.gz",reads_loc=config["reads_loc"]),
        monoReads_R2=expand("{reads_loc}/{{monos}}.2.fq.gz",reads_loc=config["reads_loc"])
    output:
        monoReads_R1=expand("{output_dir}/Preprocessing/samples/{{monos}}.1.fq.gz",output_dir=config["output_dir"]),
        monoReads_R2=expand("{output_dir}/Preprocessing/samples/{{monos}}.2.fq.gz",output_dir=config["output_dir"])
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

rule check_monos:
    input:
        monoReads_R1=expand("{output_dir}/Preprocessing/samples/{monos}.1.fq.gz",output_dir=config["output_dir"],monos=MONOS),
        monoReads_R2=expand("{output_dir}/Preprocessing/samples/{monos}.2.fq.gz",output_dir=config["output_dir"],monos=MONOS),
    output:
        checkMonos=expand("{output_dir}/Preprocessing/monoCheck.txt",output_dir=config["output_dir"])
    shell:
        """
        zcat {input.monoReads_R1} | head -n 1 > {output.checkMonos}
        """

rule check_nonMonos:
    input:
        monoReads_R1=expand("{output_dir}/Preprocessing/samples/{nonmonos}.1.fq.gz",output_dir=config["output_dir"],monos=NONMONOS),
        monoReads_R2=expand("{output_dir}/Preprocessing/samples/{nonmonos}.2.fq.gz",output_dir=config["output_dir"],monos=NONMONOS),
    output:
        checkNonMonos=expand("{output_dir}/Preprocessing/nonMonoCheck.txt",output_dir=config["output_dir"])
    shell:
        """
        zcat {input.monoReads_R1} | head -n 1 > {output.checkMonos}
        """


rule createRef:
    input:
        indRef=expand("{ref_loc}/{monos}.fa",ref_loc=config["ref_loc"],monos=MONOS),
        checkMonos=expand("{output_dir}/Preprocessing/monoCheck.txt",output_dir=config["output_dir"])
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
        