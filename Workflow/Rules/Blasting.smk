# Authors v1.0 (Legacy):    ..., ..., ... & ...
# Authors v2.0 (New):       Jelle van der Heide
# Date:                     ../../..
# 
# Considering the chance the reference is at this point still potentially contaminated with non-eukaryotic 
# DNA, the de-novo reference reads are to be blasted against the NR/NT database. The rules defined below
# execute said blast, as well as analyse it's results, filtering out unwanted reads (non-eukaryotic) accordingly.
# Besides filtering, these rules generate a number of files, which offer insight into the degree of contamination
# and what data the filtering was based upon.
# --------------------------------------------------------------------------------------------------------------------



# This rule handles blasting the reference reads against the NR/NT database. At the moment this is fully executed
# locally, which means the database could become outdated. Be sure to check with the system moderator whether the correct
# version of the database is being used. From the BLAST-results, the following columns are saved:
#       *   6 
#       *   qseqid      the query sequence id
#       *   sseqid      the target sequence id
#       *   pident      the percentage of positions that were identical
#       *   evalue      the expect value for the hit
#       *   bitscore    the bitscore for the hit
#       *   sskingdom   the super kingdom of the hit organism of origin
#       *   sscinames   the scientific name of the hit organism of origin
#       *   length      the length of the alignment
#       *   sstart      the start of the aignment within the subject
#       *   send        the end of the alignment within the subject
# -----     
# Input:    - A de-novo assembled meta-reference file.
# Output:   - The results from the BLAST query, in .tsv format
rule blast:
    params:
        blastDB=config["blastDB"]
    input:
        ref=expand("{tmp_dir}/Reference_creation/Renamed/{{sample}}.renamed.fa", tmp_dir=config["tmp_dir"])
    output:
        blastresults=temp(expand("{output_dir}/Blasting/results/{{sample}}.blastresults.tsv", output_dir=config["output_dir"]))
    log: 
        expand("{output_dir}/Logs/Blasting/{{sample}}_blast.log", output_dir=config["output_dir"])
    benchmark: 
        "../Benchmarks/blast{sample}.benchmark.tsv"
    conda: 
        "../Envs/blast.yaml"
    threads: 
        16
    resources:
        mem_mb= 500000,
        runtime= 2880,
        cpus_per_task= 16
    shell:
        """
        export BLASTDB={params.blastDB}
        blastn \
            -query {input.ref} \
            -db {params.blastDB}/nt \
            -out {output.blastresults} \
            -max_target_seqs 1 \
            -num_threads {threads} \
            -outfmt '6 qseqid sseqid pident evalue bitscore sskingdom sscinames length sstart send' \
            2> {log}
        """

rule temp_blastref:
    input:
        ref=expand("{tmp_dir}/Reference_creation/Renamed/{{sample}}.renamed.fa", tmp_dir=config["tmp_dir"]),
        blastresults=expand("{output_dir}/Blasting/results/{{sample}}.blastresults.tsv", tmp_dir=config["output_dir"])
    output:
        eukaryota_ref=expand("{output_dir}/Blasting/ref/{{sample}}/Eukaryota_ref.fa", tmp_dir=config["output_dir"])
    params:
        output_dir=config["tmp_dir"],
        sup_dir=config["sup_dir"],
        sample=expand("{{sample}}")
    conda: "../Envs/blastparse.yaml"
    resources:
        mem_mb= 10000,
        runtime= 60,
        cpus_per_task= 1    
    shell:
        """
        python Scripts/temp_solution_blast_parse.py \
            -i {input.blastresults} \
            -r {input.ref} \
            -F {params.sup_dir}/OTHER_FUNGI_NAMES.txt \
            -f {params.sup_dir}/Flowering_plant_genera_list.csv \
            -b {params.sup_dir}/Bryophytes_genera_list.csv \
            -G {params.sup_dir}/Gymnosperms_genera_list.csv \
            -P {params.sup_dir}/Pteridophytes_genera_list.csv \
            -dir '{params.output_dir}/Blasting/ref/{params.sample}' 
        """

rule blast_out:
    params:
        inputprefix=expand("{tmp_dir}/Reference_creation/Renamed/",  tmp_dir=config["tmp_dir"])
    input:
        blastresults=expand("{output_dir}/Blasting/results/{sample}.blastresults.tsv", tmp_dir=config["output_dir"],sample=MONOS)
    output:
        ref=expand("{output_dir}/Blasting/blastresults.tsv", output_dir=config["output_dir"])
    shell:
        """
        cat {params.inputprefix}*.renamed.fa > {output.ref}
        """


rule ref_out:
    params:
        inputprefix=expand("{output_dir}/Blasting",  tmp_dir=config["output_dir"])
    input:
        eukaryota_ref=expand("{output_dir}/Blasting/ref/{sample}/Eukaryota_ref.fa", tmp_dir=config["output_dir"],sample=MONOS)
    output:
        ref=expand("{output_dir}/Blasting/Eukaryota_ref.fa", output_dir=config["output_dir"])
    #log: NULL
    benchmark:
        "../Benchmarks/ref_out.benchmark.tsv"
    resources:
        mem_mb= 10000,
        runtime= 20,
        cpus_per_task= 1
    shell:
        """
        cat {params.inputprefix}/*/Eukaryota_ref.fa > {output.ref}
        """


