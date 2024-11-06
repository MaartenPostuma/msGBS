rule stats:
    input:
        bamOut=expand("{path}/mapping/mapping.bam",path=config["output_dir"])
    output:
        statscsv=expand("{path}/stats/stats.csv",path=config["output_dir"])
    shell:
        "python src/scripts/msGBS_STATS.py "
        "-i {input.bamOut} "
        "-o {output.statscsv}"
