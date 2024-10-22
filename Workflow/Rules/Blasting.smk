rule blast:
    input:
        ref=expand("{path}/output_denovo/ref.fa",path=config["output_dir"])
    output:
        blast_file=expand("{path}/output_blast/outputblast_kingdoms.tsv",path=config["output_dir"])
    params:
        blastDB=config["blastDB"]
    threads:40
    conda: "../Envs/blast.yaml"
    shell:
        """
        blastn -query {input.ref} -db {params.blastDB}/nt -out {output.blast_file} \
        -num_alignments 1 -num_threads {threads} \
        -outfmt '6 qseqid sseqid pident evalue bitscore sskingdom sscinames length sstart send'
        """


rule blastref:
    input:
        ref=expand("{path}/output_denovo/ref.fa",path=config["output_dir"]),
        blast_file=expand("{path}/output_blast/outputblast_kingdoms.tsv",path=config["output_dir"])
    output:
        refBlasted=expand("{path}/output_denovo/refBlasted.fa",path=config["output_dir"])
    params:
        outputDir=config["output_dir"]
    shell:
        "python src/scripts/blastN_parse_ref_msGBS.py  "
        "-i {input.blast_file} -r {input.ref} "
        "-F src/blastN_parse_ref_genus_lists/OTHER_FUNGI_NAMES.txt "
        "-f src/blastN_parse_ref_genus_lists/Flowering_plant_genera_list.csv "
        "-b src/blastN_parse_ref_genus_lists/Bryophytes_genera_list.csv "
        "-G src/blastN_parse_ref_genus_lists/Gymnosperms_genera_list.csv "
        "-P src/blastN_parse_ref_genus_lists/Pteridophytes_genera_list.csv "
        "-dir {params.outputDir} "


