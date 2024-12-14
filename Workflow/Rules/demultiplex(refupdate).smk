param_oligo=config["param_demultiplex"]["nOligo"]

#We always use 3 as the number of oligo for our adapters. (If we would ever change this we can account for this)
def getParam_oligo(param_oligo):
	if param_oligo == "default" or param_oligo == "":
		id = 3
	else:
		id = param_oligo
	return id

#Clone_filter removes PCR duplicates
#It compares reads that are removes reads that are completely identical 
#including the UMI/wobble/oligo's (random bases that are added at the beginning of the reads)
#as if all of these are identical it will most likely be a PCR artifact

rule polyG:
	input:
		R1=expand("{path}/{{run}}_R1.fq.gz",path=config["inputDir"]),
		R2=expand("{path}/{{run}}_R2.fq.gz",path=config["inputDir"])
	output:
		R1=temp(expand("{path}/demultiplex/trim/{{run}}_R1.fq.gz",path=config["outputDir"])),
		R2=temp(expand("{path}/demultiplex/trim/{{run}}_R2.fq.gz",path=config["outputDir"]))
	resources:
		mem_mb= 10000,
		runtime= 120,
		cpus_per_task= 6
	threads: 6
	conda:
		"env/fastp.yaml"
	shell:
		"fastp -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} --trim_poly_g"


rule clone_filter:
	input:
		barcodes=expand("{path}/{bar}", path=config["inputDir"], bar=config["barcodeFile"]),
		R1=expand("{path}/demultiplex/trim/{{run}}_R1.fq.gz",path=config["outputDir"]),
		R2=expand("{path}/demultiplex/trim/{{run}}_R2.fq.gz",path=config["outputDir"])
	params:
		outputdir=expand("{path}/demultiplex", path=config["outputDir"]),
		param_oligo=getParam_oligo(param_oligo)
	output:
		R1=expand("{path}/demultiplex/clone_filter/{{run}}_R1.1.fq.gz",path=config["outputDir"]),
		R2=expand("{path}/demultiplex/clone_filter/{{run}}_R2.2.fq.gz",path=config["outputDir"])
	conda:
		"env/stacks.yaml"
	threads: 1
	resources:
		mem_mb=lambda wc, input: max(2.5 * input.size_mb,2000),
		runtime= 6*60,
		cpus_per_task= 1,

	shell: 
		"clone_filter -1 {input.R1} -2 {input.R2} -o {params.outputdir}/clone_filter/ --oligo_len_1 {params.param_oligo} --oligo_len_2 {params.param_oligo} --inline_inline -i gzfastq"

#Stacks and the rest of the pipeline need to have specific files for barcodes, 
#the format is different and we add the control nucleotide and we need to split it per run so we can demultiplex them in parallel
#popmap (to which population do samples belong),
#and an input file for SNPFilter report to add nice colours to the plots based on a priori clustering
rule make_stacks_files:
	input:
		barcodes=expand("{path}/{bar}", path=config["inputDir"], bar=config["barcodeFile"])
	output:
		popmap=expand("{path}/stacksFiles/popmap.tsv", path=config["outputDir"]),
		popmapSNPFilter=expand("{path}/stacksFiles/SNPFilterPopMap.tsv",path=config["outputDir"]),
		barcodes=expand("{path}/stacksFiles/barcodeStacks{run}.tsv", path=config["outputDir"], bar=config["barcodeFile"],run=RUN),
		popmapSep=expand("{path}/stacksFiles/popmap{run}.tsv", path=config["outputDir"], bar=config["barcodeFile"],run=RUN)
	params:
		outputDir=expand("{path}/stacksFiles",path=config["outputDir"])
	resources:
		mem_mb= 1000,
		runtime= 1,
		cpus_per_task= 1
	conda:
		"env/R.yaml"
	shell:
		"Rscript src/demultiplex/createFilesFromBarcodeFile.R {input.barcodes} {params.outputDir}"

#We then do some snakemake magic to run process radtags for each run.
#https://stackoverflow.com/questions/4113581/snakemake-best-practice-for-demultiplexing <- based on this
# A lot of this happens in the Snakefile where a dictonary is created from the barcode file specifying in which run each sample is located.



# We than run each proces_radtags for each demultiplexing instance. The output is the process_radtags logs and a temporrary "hack directory".
# The first is a fail safe to ensure that there is actually demutliplexed data. The latter is used to generate the wild cards that are necessary for the moving/catting of the samples
# TODO the hack directory should not be necesarry but I do not know how to include something from the config in lines (84 and 102): lambda w: f"demux_tmp_{SAMPLES[w.sample]}",


rule process_radtags:
	input:
		R1=expand("{path}/demultiplex/clone_filter/{{run}}_R1.1.fq.gz",path=config["outputDir"]),
		R2=expand("{path}/demultiplex/clone_filter/{{run}}_R2.2.fq.gz",path=config["outputDir"]),
		barcodes=expand("{path}/stacksFiles/barcodeStacks{{run}}.tsv", path=config["outputDir"], bar=config["barcodeFile"])
	output:
		hackDir=directory(temp("demux_tmp_{run}")),
		log=expand("{path}/demultiplex/logs/{{run}}/process_radtags.clone_filter.log",path=config["outputDir"])
	params:
		outputDir=expand("{path}/demultiplex/logs/{{run}}/",path=config["outputDir"]),
		truncateLength=config["truncateLength"]
	conda:
		"env/stacks.yaml"
	threads: THREADSPERRUN
	resources:
		mem_mb= 10000,
		runtime= 6*60,
		cpus_per_task= 1
	shell:
		"""
		process_radtags  -1 {input.R1} -2 {input.R2} -o {params.outputDir} -b {input.barcodes} --renz_1 aseI --renz_2 nsiI -c --inline-inline --threads 1 -t {params.truncateLength}
		mkdir {output.hackDir}
		"""
#If there are no duplicates we simply move the read files from the demultiplex/logs/ directory to the demultiplex/samples folder.

#The next four rules remove individuals to which < 1000 reads were assigned to ensure the pipeline keeps running after 
#Demultiplexing
rule process_logs:
	input:
		log=expand("{path}/demultiplex/logs/{{run}}/process_radtags.clone_filter.log",path=config["outputDir"])
	output:
		log=expand("{path}/demultiplex/logs/{{run}}/perInd.tsv",path=config["outputDir"])
	conda:
		"env/stacks.yaml"
	resources:
		mem_mb= 1000,
		runtime= 1,
		cpus_per_task= 1
	shell:
		"stacks-dist-extract {input.log} per_barcode_raw_read_counts > {output.log}"

rule combine_logs:
	input:
		log=expand("{path}/demultiplex/logs/{run}/perInd.tsv",path=config["outputDir"],run=RUN)
	output:
		log=expand("{path}/demultiplex/logs/nReadsPerSample.log",path=config["outputDir"])
	shell:
		"""cat <(head -n {input.log}) <(cat {input.log} | grep -v "RAD Cutsite Not Found") > {output.log}"""

rule get_low_indvs:
	input:
		log=expand("{path}/demultiplex/logs/{{run}}/perInd.tsv",path=config["outputDir"])
	output:
		log=expand("{path}/demultiplex/logs/{{run}}/removeInds.tsv",path=config["outputDir"])
	resources:
		mem_mb= 1000,
		runtime= 1,
		cpus_per_task= 1
	shell:
		"""cat {input.log} | awk '$8<1000 {{print $2}}' > {output.log}"""

rule filter_popmap:
	input:
		popmap=expand("{path}/stacksFiles/popmap{{run}}.tsv", path=config["outputDir"]),
		removeIndv=expand("{path}/demultiplex/logs/{{run}}/removeInds.tsv",path=config["outputDir"])
	output:
		popmap=expand("{path}/stacksFiles/popmap{{run}}Filt.tsv", path=config["outputDir"])
	resources:
		mem_mb= 1000,
		runtime= 1,
		cpus_per_task= 1
	shell:
		"cat {input.popmap} | grep -f {input.removeIndv} -v	> {output.popmap}"

rule combinePerRunPopmap:
	input:
		popmap=expand("{path}/stacksFiles/popmap{run}Filt.tsv", path=config["outputDir"],run=RUN)
	output:
		popmap=expand("{path}/stacksFiles/popmapFiltDemulti.tsv", path=config["outputDir"])
	resources:
		mem_mb= 1000,
		runtime= 1,
		cpus_per_task= 1

	shell:
		"cat {input.popmap} | sort | uniq > {output.popmap}"


if DUPES==False:
	rule move_samples:
		input:
			lambda w: f"demux_tmp_{SAMPLES[w.sample]}",
		output:
			samplesR1=expand("{path}/demultiplex/samples/{{sample}}.1.fq.gz",path=config["outputDir"]),
			samplesR2=expand("{path}/demultiplex/samples/{{sample}}.2.fq.gz",path=config["outputDir"])
		params:
			outputDir=expand("{path}/demultiplex/logs/",path=config["outputDir"]),
		resources:
			mem_mb= 1000,
			runtime= 1,
			cpus_per_task= 1

		shell:
			"""
			mv {params.outputDir}/*/{wildcards.sample}.1.fq.gz {output.samplesR1}
			mv {params.outputDir}/*/{wildcards.sample}.2.fq.gz {output.samplesR2}
			"""
#If there are duplicates in the log file we cannot move the files. As two files with the same name in different runs should be catted
#For whatever reason snakemake selects only one combination of run/sample from the dictionary if there are two or more? So by including the logs files as inputs we ensure all
#instances of process_radtags are actually done running.
#When they are done all samples with the same name are catted together and the files in the logs directory are removed.
if DUPES==True:
	rule cat_samples:
		input:
			lambda w: f"demux_tmp_{SAMPLES[w.sample]}",
			log=expand("{path}/demultiplex/logs/{run}/process_radtags.clone_filter.log",path=config["outputDir"],run=RUN)
		output:
			samplesR1=expand("{path}/demultiplex/samples/{{sample}}.1.fq.gz",path=config["outputDir"]),
			samplesR2=expand("{path}/demultiplex/samples/{{sample}}.2.fq.gz",path=config["outputDir"])
		params:
			outputDir=expand("{path}/demultiplex/logs/",path=config["outputDir"])
		resources:
			mem_mb= 1000,
			runtime= 1,
			cpus_per_task= 1
		shell:
			"""
			cat {params.outputDir}*/{wildcards.sample}.1.fq.gz > {output.samplesR1}
			cat {params.outputDir}*/{wildcards.sample}.2.fq.gz > {output.samplesR2}
			rm {params.outputDir}*/{wildcards.sample}.1.fq.gz
			rm {params.outputDir}*/{wildcards.sample}.2.fq.gz
			rm {params.outputDir}*/{wildcards.sample}.rem.1.fq.gz
			rm {params.outputDir}*/{wildcards.sample}.rem.2.fq.gz
			"""
