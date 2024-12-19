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
        ref=expand("{ref_dir}/ref.fa", ref_dir=config["ref_dir"])
    output:
        blastresults=expand("{blasting_out}/blastresults.tsv", blasting_out=config["blasting_out"])
    log: 
        "../Logs/Blasting/blast.log"
    benchmark: 
        "../Benchmarks/blast.benchmark.tsv"
    conda: 
        "../Envs/blast.yaml"
    threads: 
        32
    shell:
        """
        export BLASTDB={params.blastDB}
        blastn \
            -query {input.ref} \
            -db {params.blastDB}/nt \
            -out {output.blastresults} \
            -num_alignments 1 \
            -num_threads {threads} \
            -outfmt '6 qseqid sseqid pident evalue bitscore sskingdom sscinames length sstart send' \
            2> {log}
        """

rule temp_blastref:
    input:
        ref=expand("{ref_dir}/ref.fa", ref_dir=config["ref_dir"]),
        blastresults=expand("{blasting_out}/blastresults.tsv", blasting_out=config["blasting_out"])
    output:
        eukaryota_ref=expand("{blasting_out}/Eukaryota_ref.fa", blasting_out=config["blasting_out"]),
    params:
        blasting_out=config["blasting_out"],
        sup_dir=config["sup_dir"]
    conda: "../Envs/blastparse.yaml"
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
            -dir '{params.blasting_out}' 
        """

"""
# While this rule may seem a bit large at first, what it does really isn't all that impressive. In short, this rule parses the
# BLAST-result file and extracts data of interest. The main focus here lies on filtering non-eukaryotic reads from the reference.
# Reads that belong to other kingdoms are saved in separate ref files, as well as their scientific names and headers (both also separate). 
# Lastly, an overview is generated for the number of reads that have been assigned to which super kingdom.
# -----     
# Input:    - The de-novo assembled meta-reference read-file.
#           - The BLAST-result file.
#           - The list of known genusses, which has potentially been generated in a prior run.
# Output:   - Files that contain the names found within the hits of the BLAST results that belong to a specific super kingdom.
#           - Files that contain the headers of above mentioned hits for a specific super kingdom.
#           - Files that contain the filtered meta-reference for a specific super kingdom.
#           - The statistics of the blast query, specifically how many reads were assigned to which super kingdom.
rule blastref:
    #params NULL
    input:
        ref=expand("{ref_dir}/ref.fa", ref_dir=config["ref_dir"]),
        blastresults=expand("{blasting_out}/blastresults.tsv", blasting_out=config["blasting_out"]),
        genus_lists=expand("{sup_dir}/{genus}", sup_dir=config["sup_dir"], genus=config["genus"])
    output:
        bryophyta_names=expand("{blasting_out}/Bryophyta_names.txt", blasting_out=config["blasting_out"]),
        eukaryota_names=expand("{blasting_out}/Eukaryota_names.txt", blasting_out=config["blasting_out"]),
        imperia_names=expand("{blasting_out}/Imperia_names.txt", blasting_out=config["blasting_out"]),
        fungi_names=expand("{blasting_out}/Fungi_names.txt", blasting_out=config["blasting_out"]),
        arbuscular_names=expand("{blasting_out}/Arbuscular_names.txt", blasting_out=config["blasting_out"]),
        gymnospermae_names=expand("{blasting_out}/Gymnospermae_names.txt", blasting_out=config["blasting_out"]),
        pteridophyta_names=expand("{blasting_out}/Pteridophyta_names.txt", blasting_out=config["blasting_out"]),
        angiospermae_names=expand("{blasting_out}/Angiospermae_names.txt", blasting_out=config["blasting_out"]),
        archaea_names=expand("{blasting_out}/Archaea_names.txt", blasting_out=config["blasting_out"]),
        bacteria_names=expand("{blasting_out}/Bacteria_names.txt", blasting_out=config["blasting_out"]),
        bryophyta=expand("{blasting_out}/Bryophyta.txt", blasting_out=config["blasting_out"]),
        eukaryota=expand("{blasting_out}/Eukaryota.txt", blasting_out=config["blasting_out"]),
        imperia=expand("{blasting_out}/Imperia.txt", blasting_out=config["blasting_out"]),
        fungi=expand("{blasting_out}/Fungi.txt", blasting_out=config["blasting_out"]),
        arbuscular=expand("{blasting_out}/Arbuscular.txt", blasting_out=config["blasting_out"]),
        gymnospermae=expand("{blasting_out}/Gymnospermae.txt", blasting_out=config["blasting_out"]),
        pteridophyta=expand("{blasting_out}/Pteridophyta.txt", blasting_out=config["blasting_out"]),
        angiospermae=expand("{blasting_out}/Angiospermae.txt", blasting_out=config["blasting_out"]),
        archaea=expand("{blasting_out}/Archaea.txt", blasting_out=config["blasting_out"]),
        bacteria=expand("{blasting_out}/Bacteria.txt", blasting_out=config["blasting_out"]),
        eukaryota_ref=expand("{blasting_out}/Eukaryota_ref.fa", blasting_out=config["blasting_out"]),
        imperia_ref=expand("{blasting_out}/Imperia_ref.fa", blasting_out=config["blasting_out"]),
        fungi_ref=expand("{blasting_out}/Fungi_ref.fa", blasting_out=config["blasting_out"]),
        arbuscular_ref=expand("{blasting_out}/Arbuscular_ref.fa", blasting_out=config["blasting_out"]),
        angiospermae_ref=expand("{blasting_out}/Angiospermae_ref.fa", blasting_out=config["blasting_out"]),
        archaea_ref=expand("{blasting_out}/Archaea_ref.fa", blasting_out=config["blasting_out"]),
        bacteria_ref=expand("{blasting_out}/Bacteria_ref.fa", blasting_out=config["blasting_out"]),
        stats=expand("{blasting_out}/blastStatistics.txt", blasting_out=config["blasting_out"])
    #log: NULL
    benchmark:
        "../Benchmarks/blastref.benchmark.tsv"
    #conda NULL
    threads: 
        32
    shell:      
        """"""
        touch {output.bryophyta_names}
        touch {output.eukaryota_names}
        touch {output.imperia_names}
        touch {output.fungi_names}
        touch {output.arbuscular_names}
        touch {output.gymnospermae_names}
        touch {output.pteridophyta_names}
        touch {output.angiospermae_names}
        touch {output.archaea_names}
        touch {output.bacteria_names}
        touch {output.bryophyta}
        touch {output.eukaryota}
        touch {output.imperia}
        touch {output.fungi}
        touch {output.arbuscular}
        touch {output.gymnospermae}
        touch {output.pteridophyta}
        touch {output.angiospermae}
        touch {output.archaea}
        touch {output.bacteria}
        touch {output.eukaryota_ref}
        touch {output.imperia_ref}
        touch {output.fungi_ref}
        touch {output.arbuscular_ref}
        touch {output.angiospermae_ref}
        touch {output.archaea_ref}
        touch {output.bacteria_ref}

        Bryophyta=()            #Mosses and worts
        Eukaryota=()            #Eukaryota
        Imperia=()              #Viruses
        Fungi=()                #Fungi
        Arbuscular=()           #Arbuscular mycorrhiza (soil borne fungi)
        Gymnospermae=()         #Seed producing  (unenclosed seeds)
        Pteridophyta=()         #Spore producing vascular plants
        Angiospermae=()         #Flowering plants
        Archaea=()              #Archaea
        Bacteria=()             #Bacteria

        Bryophyta_contigs=()            
        Eukaryota_contigs=()            
        Imperia_contigs=()              
        Fungi_contigs=()           
        Arbuscular_contigs=()
        Gymnospermae_contigs=()
        Pteridophyta_contigs=()
        Angiospermae_contigs=()
        Archaea_contigs=()
        Bacteria_contigs=()

        echo "Total number of blast hits: " > {output.stats}
        awk 'END {{ print NR}}' {input.ref} > {output.stats}
        echo "-----------------------------------------" >> {output.stats}
        echo "Total number of blast hits per Kingdom/phylum of interest:" >> {output.stats}

        while IFS=";" read -r Br Eu Im Fu Gy Pt An Am; do
            if [[ $Br = *[![:space:]]* ]];
            then
                Bryophyta+=("${{Br//}}")
            fi
            if [[ $Eu = *[![:space:]]* ]];
            then
                Eukaryota+=("${{Eu//}}")
            fi
            if [[ $Im = *[![:space:]]* ]];
            then
                Imperia+=("${{Im//}}")
            fi
            if [[ $Fu = *[![:space:]]* ]];
            then
                Fungi+=("${{Fu//}}")
            fi
            if [[ $Gy = *[![:space:]]* ]];
            then
                Gymnospermae+=("${{Gy//}}")
            fi
            if [[ Pt = *[![:space:]]* ]];
            then
                Pteridophyta+=("${{Pt//}}")
            fi
            if [[ $An = *[![:space:]]* ]];
            then
                Angiospermae+=("${{An//}}")
            fi
            if [[ $Am = *[![:space:]]* ]];
            then
                Arbuscular+=("${{Am//}}")
            fi
        done < {input.genus_lists}
    
        while IFS='	' read -r contig hit pident e_value bp_hit kingdoms fullname alignment_length start stop; do
         e_value_min=${{e_value#*-}}    
         add=false
         if [[ $e_value_min != *"."* ]]; then
             if [ $alignment_length -gt 40 ] && [ $e_value_min -gt 20 ]; then
                 add=true
             fi 
         fi
        case $kingdoms in
            Bacteria)  
                if [ "$add" = true ]; then
                    Bacteria_contigs+=("$contig")
                    echo "$fullname" >>  {output.bacteria_names} 
                fi
                ;;
            Archaea)
                if [ "$add" = true ]; then
                    Archaea_contigs+=("$contig")
                    echo "$fullname" >>  {output.archaea_names} 
                fi
                ;;    
            Viruses)
                if [ "$add" = true ]; then
                    Imperia_contigs+=("$contig")
                    echo "$fullname" >>  {output.imperia_names}
                fi
                ;;    
            N/A)
                if [ "$add" = true ]; then
                    Eukaryota_contigs+=("$contig")
                    echo "$fullname" >>  {output.eukaryota_names}
                fi
                ;;
            Eukaryota)
                genera=($fullname)
                main_genera="${{genera[@]:0:1}}"
                sub_genera="${{genera[@]:0:2}}" 
                if [[ ${{Eukaryota[@]}} =~ $main_genera ]] && [ "$add" = true ]; then
                    Eukaryota_contigs+=("$contig")
                    echo "$fullname" >>  {output.eukaryota_names}
                fi
                if [[ ${{Arbuscular[@]}} =~ $main_genera ]] && [ "$add" = true ]; then
                    Arbuscular_contigs+=("$contig")
                    Fungi_contigs+=("$contig")
                    echo "$sub_genera" >> {output.arbuscular_names}
                    echo "$sub_genera" >> {output.fungi_names} 
                elif [[ ${{Fungi[@]}} =~ $main_genera ]] && [ "$add" = true ]; then
                    Fungi_contigs+=("$contig")
                    echo "$sub_genera" >> {output.fungi_names}
                elif [[ ${{Angiospermae[@]}} =~ $main_genera ]] && [ "$add" = true ]; then
                    Angiospermae_contigs+=("$contig")
                    echo "$fullname" >> {output.angiospermae_names}
                elif [[ ${{Bryophyta[@]}} =~ $main_genera ]] && [ "$add" = true ]; then
                    Bryophyta_contigs+=("$contig")
                    echo "$fullname" >> {output.bryophyta_names}
                elif [[ ${{Gymnospermae[@]}} =~ $main_genera ]] && [ "$add" = true ]; then
                    Gymnospermae_contigs+=("$contig")
                    echo "$fullname" >> {output.gymnospermae_names}
                elif [[ ${{Pteridophyta[@]}} =~ $main_genera ]] && [ "$add" = true ]; then
                    Pteridophyta_contigs+=("$contig")
                    echo "$fullname" >> {output.pteridophyta_names}
                fi
                ;;
            *) 
                if [ "$add" = true ]; then
                    Eukaryota_contigs+=("$contig")
                    echo "$fullname" >>  {output.eukaryota_names}
                fi
                ;;
        esac
        done < {input.blastresults}

        echo "\tBryophyta:\t\t" "${{#Bryophyta_contigs[@]}}" >> {output.stats}
        echo "\tEukaryota:\t\t" "${{#Eukaryota_contigs[@]}}"  >> {output.stats}
        echo "\tImperia:\t\t" "${{#Imperia_contigs[@]}}"  >> {output.stats}
        echo "\tFungi:\t\t\t" "${{#Fungi_contigs[@]}}"  >> {output.stats}
        echo "\tArbuscular:\t\t" "${{#Arbuscular_contigs[@]}}"  >> {output.stats}
        echo "\tGymnospermae:\t" "${{#Gymnospermae_contigs[@]}}"  >> {output.stats}
        echo "\tPteridophyta:\t" "${{#Pteridophyta_contigs[@]}}"  >> {output.stats}
        echo "\tAngiospermae:\t" "${{#Angiospermae_contigs[@]}}"  >> {output.stats}
        echo "\tArchaea:\t\t" "${{#Archaea_contigs[@]}}"  >> {output.stats}
        echo "\tBacteria:\t\t" "${{#Bacteria_contigs[@]}}" >> {output.stats}                          

        for line in ${{Bryophyta_contigs[@]}}
        do
            echo "$line" >> {output.bryophyta}
        done
        for line in ${{Eukaryota_contigs[@]}}
        do
            echo "$line" >> {output.eukaryota}
        done
        for line in ${{Imperia_contigs[@]}}
        do
            echo "$line" >> {output.imperia}
        done
        for line in ${{Fungi_contigs[@]}}
        do
            echo "$line" >> {output.fungi}
        done
        for line in ${{Arbuscular_contigs[@]}}
        do
            echo "$line" >> {output.arbuscular}
        done
        for line in ${{Gymnospermae_contigs[@]}}
        do
            echo "$line" >> {output.gymnospermae}
        done
        for line in ${{Pteridophyta_contigs[@]}}
        do
            echo "$line" >> {output.pteridophyta}
        done
        for line in ${{Angiospermae_contigs[@]}}
        do
            echo "$line" >> {output.angiospermae}
        done  
        for line in ${{Archaea_contigs[@]}}
        do
            echo "$line" >> {output.archaea}
        done
        for line in ${{Bacteria_contigs[@]}}
        do
            echo "$line" >> {output.bacteria}
        done        

        echo "-----------------------------------------" >> {output.stats}
        echo "Total number of lines in reference-sequence before blast-filter:" >> {output.stats}
        awk 'END {{ print NR}}' {input.ref} >> {output.stats}
        switchcase="NONE"
        while read line; do
            if [[ $line == \>* ]]; then
                #line="${{line:1}}"    # deze line moet er uit om bowtie te kunnen runnen
                if [[ ${{Angiospermae_contigs[@]}} =~ $line ]]; then
                    echo "$line" >>  {output.angiospermae_ref}
                    Angiospermae_contigs=( "${{Angiospermae_contigs[@]/$line}}" )
                    echo "$line" >>  {output.eukaryota_ref} 
                    Eukaryota_contigs=( "${{Eukaryota_contigs[@]/$line}}" )
                    switchcase="ANGIOSPERMAE"
                elif [[ ${{Eukaryota_contigs[@]}} =~ $line ]]; then
                    echo "$line" >>  {output.eukaryota_ref}
                    Eukaryota_contigs=( "${{Eukaryota_contigs[@]/$line}}" )
                    switchcase="EUKARYOTA"
                elif [[ ${{Bacteria_contigs[@]}} =~ $line ]]; then
                    echo "$line" >>  {output.bacteria_ref} 
                    Bacteria_contigs=( "${{Bacteria_contigs[@]/$line}}" )
                    switchcase="BACTERIA"
                elif [[ ${{Arbuscular_contigs[@]}} =~ $line ]]; then
                    echo "$line" >>  {output.arbuscular_ref}
                    Arbuscular_contigs=( "${{Arbuscular_contigs[@]/$line}}" )
                    switchcase="ARBUSCULAR"
                elif [[ ${{Fungi_contigs[@]}} =~ $line ]]; then
                    echo "$line" >>  {output.fungi_ref}
                    Fungi_contigs=( "${{Fungi_contigs[@]/$line}}" )
                    switchcase="FUNGI"
                elif [[ ${{Archaea_contigs[@]}} =~ $line ]]; then
                    echo "$line" >>  {output.archaea_ref}
                    Archaea_contigs=( "${{Archaea_contigs[@]/$line}}" )
                    switchcase="ARCHAEA"
                elif [[ ${{Imperia_contigs[@]}} =~ $line ]]; then
                    echo "$line" >>  {output.imperia_ref} 
                    Imperia_contigs=( "${{Imperia_contigs[@]/$line}}" )
                    switchcase="IMPERIA"
                else
                    echo "$line" >>  {output.eukaryota_ref}
                    Eukaryota_contigs=( "${{Eukaryota_contigs[@]/$line}}" )
                    switchcase="EUKARYOTA"
                fi
            else 
                case $switchcase in
                    ANGIOSPERMAE)
                        echo "$line" >>  {output.angiospermae_ref}
                        echo "$line" >>  {output.eukaryota_ref}  
                        ;;
                    EUKARYOTA)
                        echo "$line" >>  {output.eukaryota_ref}
                        ;;
                    BACTERIA) 
                        echo "$line" >>  {output.bacteria_ref} 
                        ;;
                    ARBUSCULAR)
                        echo "$line" >>  {output.arbuscular_ref}
                        ;;
                    FUNGI)
                        echo "$line" >>  {output.fungi_ref} 
                        ;;
                    ARCHAEA)
                        echo "$line" >>  {output.archaea_ref}
                        ;;
                    IMPERIA)
                        echo "$line" >>  {output.imperia_ref} 
                        ;;
                    *)
                        echo "$line" >>  {output.eukaryota_ref}
                        ;;
                esac
            fi
        done < {input.ref}
        """