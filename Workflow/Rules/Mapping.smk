rule mapping_Bowtie2_index:
    input: 
        refBlasted=expand("{path}/output_blast/Eukaryota_ref.fa",path=config["output_dir"]),
    output:
        index=expand("{path}/mapping/index.1.bt2", path=config["output_dir"])
    conda: "../Envs/bowtie2.yaml"
    threads: 4
    shell:
        """
        bowtie2-build -f {input.refBlasted} ../Output/mapping/index -p 4
        """

rule mapping_Bowtie2:   
    params:
        sample='{sample}',
    input:
        index=expand("{path}/mapping/index.1.bt2",  path=config["output_dir"]),
        r1=expand("{tmp}/Preprocessing/Sampleheaders/{{sample}}.1.fq.gz",  tmp=config["tmp_dir"]),#, sample=SAMPLES),
        r2=expand("{tmp}/Preprocessing/Sampleheaders/{{sample}}.2.fq.gz",  tmp=config["tmp_dir"])#, sample=SAMPLES)
    conda: "../Envs/bowtie2.yaml"
    benchmark:"../Benchmarks/Bowtie2-{sample}.benchmark.tsv"
    output:
        samOut=temp(expand("{path}/mapping/mapping_sq_{{sample}}.sam",path=config["output_dir"])),
        bamOut=expand("{path}/mapping/mapping_sq_{{sample}}.bam",path=config["output_dir"])
    threads: 4
    shell:
        """
        bowtie2 -x ../Output/mapping/index -1 {input.r1} -2 {input.r2} -q --end-to-end --very-fast --threads 4 -S {output.samOut}
        samtools view -b -o {output.bamOut} {output.samOut}
        """ #k10

#rule mapping_bwa_index:
#    input: 
#        refBlasted=expand("{path}/output_blast/Eukaryota_ref.fa",path=config["output_dir"]),
#    output:
#        index=expand("{path}/mapping/index.amb", path=config["output_dir"])
#    conda: "../Envs/bwa.yaml"
#    threads: 4
#    shell:
#        """
#        bwa index -p ../Output/mapping/index {input.refBlasted}
#        touch tmpfile.txt
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
#    output:
#        samOut=temp(expand("{path}/mapping/mapping_rg_{{sample}}.sam",path=config["output_dir"])),
#        bamOut=temp(expand("{path}/mapping/mapping_rg_{{sample}}.bam",path=config["output_dir"])),
#    threads: 4
#    shell:
#        """
#        bwa mem -t 4 -c 10 -R '@RG\\tID:{params.sample}\\tSM:{params.sample}' ../Output/mapping/index {input.r1} {input.r2} > {output.samOut}
#        samtools view -b -o {output.bamOut} {output.samOut}
#        """ #c 10

#rule mapping_GEM3:

#rule mapping_Segemehl:

#rule mapping_star_index:
#    params:
#        tempMaps=expand("../Misc/mapping/indexed")
#    input:
#        refBlasted=expand("{path}/output_blast/Eukaryota_ref.fa",path=config["output_dir"]),
#    output:
#        genome=expand("{output}/mapping/index/Genome" , output=config["output_dir"]),
#    conda: "../Envs/star.yaml"
#    shell:
#        """
#        #logLen=awk'{{print log(*|||*genoomlengte*|||*)/log(2)}}'
#        #minimizedLog=loglen / 2 - 1
#        #if[minimizedLog > 14] then
#        #    genomeSAindexNbases=14
#        #else
#        #    genomeSAindexNbases= minimizedlog (splits op . en dan eerste waarde, of gewoon alleen flooren ofzo)#
#        STAR --genomeSAindexNbases 10 --runThreadN 4 --runMode genomeGenerate --genomeDir ../Output/mapping/index --genomeFastaFiles {input.refBlasted} --outTmpDir {params.tempMaps}
#        """

#rule mapping_Star:
#    params:
#        r1out=temp(expand("{tmp}/Preprocessing/Sampleheaders/{{sample}}.1.fq", tmp=config["tmp_dir"])),
#        r2out=temp(expand("{tmp}/Preprocessing/Sampleheaders/{{sample}}.2.fq", tmp=config["tmp_dir"])),
#        sample='{sample}'
#    input:
#        genome=expand("{output}/mapping/index/Genome" , output=config["output_dir"]),
#        refBlasted=expand("{path}/output_blast/Eukaryota_ref.fa",path=config["output_dir"]),
#        r1=expand("{tmp}/Preprocessing/Sampleheaders/{{sample}}.1.fq.gz",  tmp=config["tmp_dir"]),#, sample=SAMPLES),
#        r2=expand("{tmp}/Preprocessing/Sampleheaders/{{sample}}.2.fq.gz",  tmp=config["tmp_dir"])#, sample=SAMPLES)
#    conda: "../Envs/star.yaml"
#    output:        
#        #r1out=temp(expand("{tmp}/Preprocessing/Sampleheaders/{{sample}}.1.fq", tmp=config["tmp_dir"])),
#        #r2out=temp(expand("{tmp}/Preprocessing/Sampleheaders/{{sample}}.2.fq", tmp=config["tmp_dir"])),
#        samOut=temp(expand("{path}/mapping/mapping_rg_{{sample}}.sam",path=config["output_dir"])),
#        bamOut=expand("{path}/mapping/mapping_sq_{{sample}}.bam",path=config["output_dir"])
#    threads: 16
#    shell:
#        """
#	gunzip {input.r1}
#	gunzip {input.r2}
#        STAR runThreadN 16 --genomeDir ../Output/mapping/index --readFilesIn {params.r1out} {params.r2out} --outSAMattributes NM MD AS --outSAMtype SAM --outFileNamePrefix ../Output/mapping/{params.sample}_ --outFilterMatchNminOverLread 0.95 --clip3pNbases 1 1 --outSAMorder PairedKeepInputOrder --outFilterMultimapScoreRange 0 --alignEndsType Extend5pOfRead1 --scoreGapNoncan 0 --scoreGapGCAG 0 --scoreGapATAC 0 --scoreDelOpen 0 --scoreDelBase 0 --scoreInsOpen 0 --scoreInsBase 0 --alignMatesGapMax 20 --readMapNumber -1
#        mv ../Output/mapping/{params.sample}_Aligned.out.sam {output.samOut}
#        samtools view -b -o {output.bamOut} {output.samOut}
#        """

        #"python Scripts/Map_STAR_snake.py --tmpdir ../Misc/Star "
        #"--input_dir ../Misc/Preprocessing/Sampleheaders "
        #"--output_dir ../Output/mapping "
        #"--threads 4 "
        #"--barcodes ../Input/"

        #genomeSAindexNbases = min(14, math.log(genome_len, 2) / 2 - 1)
        #index_cmd += ' --genomeSAindexNbases %i' % genomeSAindexNbases
        #genomeChrBinNbits = min(18, math.log(genome_len / no_clusters, 2))
        #index_cmd += ' --genomeChrBinNbits %i' % genomeChrBinNbits
#def addRG(in_files, args):
#    """make header for output bamfile and split in watson and crick"""
    # define readgroup header lines by combining the following

#    with open(args.barcodes, 'r') as barcodes:
#        sam_out = open(in_files['header'], 'a')
#        header = barcodes.readline().split('\t')
#        for line in barcodes:
#            RG = ['@RG']
#            split_line = line.split('\t')
#            if args.species and 'Species' in header:
#                if split_line[(header.index('Species'))] != args.species:
#                    continue
#            fc = split_line[(header.index('Flowcell'))]
#            lane = split_line[(header.index('Lane'))]
#            sample = split_line[(header.index('Sample'))]
#            RG.append('ID:%s_%s_%s' % (fc, lane, sample))
#            RG.append('SM:%s' % (sample))
#            RG.append('LB:%s_%s' % (fc, sample))
#            RG.append('PL:ILLUMINA\n')
#            sam_out.write('\t'.join(RG))
#    sam_out.close()
#    return in_files


#def make_header(args):
#    """Make header for watson and crick bam file"""
#    header = os.path.join(args.output_dir, 'header.sam')
#    args.header = header
#    header_handle = open(header, 'w')
#    header_handle.write('@HD\tVN:1.4\n')
#    joined_sam = open(os.path.join(args.output_dir, 'joinedAligned.out.sam'))
#    merged_sam = open(os.path.join(args.output_dir, 'mergedAligned.out.sam'))
#    for line in joined_sam:
#        if line.startswith('@'):
#            if line.startswith('@SQ'):
#                header_handle.write(line)
#        else:
#            break
#    for line in merged_sam:
#        if line.startswith('@'):
#            if line.startswith('@SQ'):
#                header_handle.write(line)
#            elif not line.startswith('@HD'):
#                header_handle.write(line)
#        else:
#            break
#    header_handle.close()
#    in_files = {'header': os.path.join(args.output_dir, 'header.sam')}
#    addRG(in_files, args)
#    return args


#def bam_output(args):
#    """Generate watson and crick output bam file"""
#
#    merged_sam = os.path.join(args.output_dir, 'mergedAligned.out.sam')
#    joined_sam = os.path.join(args.output_dir, 'joinedAligned.out.sam')#
#
#    out_sam = tempfile.NamedTemporaryFile(prefix='tmp', suffix='.sam', dir=args.output_dir)
#    # rewrite sam file merged and joined for watson and crick
#    parse_sam(merged_sam, out_sam.name, 'merged', 'tmp')
#    # TODO: determine why joined reads have more soft-clips or single read matches
#    parse_sam(joined_sam, out_sam.name, 'joined', 'tmp')
#    # convert to sorted and indexed bam

#    cmd = 'cat %s %s |samtools view -@ 4 -Shb |sambamba sort --tmpdir %s -m 4GB -t %s -o %s  /dev/stdin' % (args.header,
#                                                                                                            out_sam.name,args.tmpdir,
#                                                                                                            args.threads,
#                                                                                                            os.path.join(
#                                                                                                                args.output_dir,
#                                                                                                                'out.bam'))
#if not os.path.exists(os.path.join(args.output_dir,'header.sam')): #main
#    args = process_reads_merged(args)
#    args = process_reads_joined(args)
#    args = index_STAR(args)
#    args = map_STAR(args)
#    args = make_header(args)
#    args = bam_output(args)
