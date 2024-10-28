#rule blast:
#    input:
#        ref=expand("{path}/output_denovo/ref.fa",path=config["output_dir"])
#    output:
#        blast_file=expand("{path}/output_blast/outputblast_kingdoms.tsv",path=config["output_dir"])
#    params:
#        blastDB=config["blastDB"]
#    threads:40
#    conda: "../Envs/blast.yaml"
#    shell:
#        """
#        blastn -query {input.ref} -db {params.blastDB}/nt -out {output.blast_file} \
#        -num_alignments 1 -num_threads {threads} \
#        -outfmt '6 qseqid sseqid pident evalue bitscore sskingdom sscinames length sstart send'
#        """

rule blastref:
    input:
        ref=expand("{path}/output_denovo/ref.fa",path=config["output_dir"]),
        blast_file=expand("{path}/output_blast/outputblast_kingdoms.tsv",path=config["output_dir"]),
        genus_lists=expand("{path}", path=config["genus"])
    output:
        refBlasted=expand("{path}/output_denovo/Eukaryota_ref.fa",path=config["output_dir"])
    params:
        outputDir=expand("{path}/output_blast/", path=config["output_dir"])
    shell:
        """
        add_check () {{
            e_value_min=${{1#*-}}
            if [[ $e_value_min != *"."* ]]; then
                if [ $2 -gt 40 ] && [ $e_value_min -gt 20 ]; then
                    return 1
                fi 
            fi
            return 0
        }}

        declare -A genus_counts=( ["Bryophyta"]=0 ["Eukaryota"]=0 ["Imperia"]=0 ["Fungi"]=0 ["Gymnospermae"]=0 ["Pteridophyta"]=0 ["Angiospermae"]=0 ["Bacteria"]=0 ["Archaea"]=0  ["N/A"]=0)

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

        echo "nr of blast hits: "
        wc -l {input.blast_file}
        echo "-----------------------------------------"

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
        case $kingdoms in
            Bacteria)
                add_check $e_value $alignment_length       
                if [ $? -eq 1 ]; then
                    Bacteria_contigs+=("$contig")
                    echo "$fullname" >>  {params.outputDir}Bacteria_names.txt  
                    ((Bacteria_count+=1))
                fi
                ;;
            Archaea)
                add_check $e_value $alignment_length       
                if [ $? -eq 1 ]; then
                    Archaea_contigs+=("$contig")
                    echo "$fullname" >>  {params.outputDir}Archaea_names.txt  
                    ((Archaea_count+=1))
                fi
                ;;    
            Viruses)
                add_check $e_value $alignment_length       
                if [ $? -eq 1 ]; then
                    Imperia_contigs+=("$contig")
                    echo "$fullname" >>  {params.outputDir}Imperia_names.txt  
                    ((Imperia_count+=1))
                fi
                ;;    
            N/A)
                add_check $e_value $alignment_length       
                if [ $? -eq 1 ]; then
                    Eukaryota_contigs+=("$contig")
                    echo "$fullname" >>  {params.outputDir}Eukaryota_names.txt  
                    ((Eukaryota_count+=1))
                fi
                ;;
            Eukaryota)
                genera=($fullname)
                main_genera="${{genera[@]:0:1}}"
                sub_genera="${{genera[@]:0:2}}"
                add_check $e_value $alignment_length   
                add=$?
                if [[ ${{Eukaryota[@]}} =~ $main_genera ]] && [ $add -eq 1 ]; then
                    Eukaryota_contigs+=("$contig")
                    echo "$fullname" >>  {params.outputDir}Eukaryota_names.txt  
                    (( genus_counts["Eukaryota"]++ ))
                fi
                if [[ ${{Arbuscular[@]}} =~ $main_genera ]] && [ $add -eq 1 ]; then
                    Arbuscular_contigs+=("$contig")
                    echo "$sub_genera" >> {params.outputDir}Arbuscular_names.txt  
                    (( genus_counts["Fungi"]++ ))
                    (( genus_counts["Arbuscular"]++ ))
                elif [[ ${{Fungi[@]}} =~ $main_genera ]] && [ $add -eq 1 ]; then
                    Fungi_contigs+=("$contig")
                    echo "$sub_genera" >> {params.outputDir}Fungi_names.txt  
                    (( genus_counts["Fungi"]++ ))
                elif [[ ${{Angiospermae[@]}} =~ $main_genera ]] && [ $add -eq 1 ]; then
                    Angiospermae_contigs+=("$contig")
                    echo "$fullname" >> {params.outputDir}Angiospermae_names.txt  
                    (( genus_counts["Angiospermae"]++ ))
                elif [[ ${{Bryophyta[@]}} =~ $main_genera ]] && [ $add -eq 1 ]; then
                    Bryophyta_contigs+=("$contig")
                    echo "$fullname" >> {params.outputDir}Bryophyta_names.txt 
                    (( genus_counts["Bryophyta"]++ ))
                elif [[ ${{Gymnospermae[@]}} =~ $main_genera ]] && [ $add -eq 1 ]; then
                    Gymnospermae_contigs+=("$contig")
                    echo "$fullname" >> {params.outputDir}Gymnospermae_names.txt 
                    (( genus_counts["Gymnospermae"]++ ))
                elif [[ ${{Pteridophyta[@]}} =~ $main_genera ]] && [ $add -eq 1 ]; then
                    Pteridophyta_contigs+=("$contig")
                    echo "$fullname" >> {params.outputDir}Pteridophyta_names.txt 
                    (( genus_counts["Pteridophyta"]++ ))
                fi
                ;;
            *) 
                add_check $e_value $alignment_length       
                if [ $? -eq 1 ]; then
                    Eukaryota_contigs+=("$contig")
                    echo "$fullname" >>  {params.outputDir}Eukaryota_names.txt  
                    (( genus_counts["Eukaryota"]++ ))
                fi
                ;;
        esac
        done < {input.blast_file}
    
        echo "perensap"


        for i in "${{!genus_counts[@]}}"
        do
            echo "$i" "${{genus_counts[$i]}}"
        done

        for line in ${{Bryophyta_contigs[@]}}
        do
            echo "$line" >> {params.outputDir}Bryophyta.txt
        done
        for line in ${{Eukaryota_contigs[@]}}
        do
            echo "$line" >> {params.outputDir}Eukaryota.txt
        done
        for line in ${{Imperia_contigs[@]}}
        do
            echo "$line" >> {params.outputDir}Imperia.txt
        done
        for line in ${{Fungi_contigs[@]}}
        do
            echo "$line" >> {params.outputDir}Fungi.txt
        done
        for line in ${{Arbuscular_contigs[@]}}
        do
            echo "$line" >> {params.outputDir}Arbuscular.txt
        done
        for line in ${{Gymnospermae_contigs[@]}}
        do
            echo "$line" >> {params.outputDir}Gymnospermae.txt
        done
        for line in ${{Pteridophyta_contigs[@]}}
        do
            echo "$line" >> {params.outputDir}Pteridophyta.txt
        done
        for line in ${{Angiospermae_contigs[@]}}
        do
            echo "$line" >> {params.outputDir}Angiospermae.txt
        done  
        for line in ${{Archaea_contigs[@]}}
        do
            echo "$line" >> {params.outputDir}Archaea.txt
        done
        for line in ${{Bacteria_contigs[@]}}
        do
            echo "$line" >> {params.outputDir}Bacteria.txt
        done        

        echo "nr of ref-lines: "
        wc -l {input.ref}
        echo "-----------------------------------------"
        switchcase="NONE"
        while read line; do
            if [[ $line == \>* ]]; then
                line="${{line:1}}"    # eigenlijk rstrip, maar hoeft denk ik niet
                if [[ ${{Angiospermae_contigs[@]}} =~ $line ]]; then
                    echo "zalmkoeskoes"
                    echo "$line" >>  {params.outputDir}Angiospermae_ref.fa
                    Angiospermae_contigs=( "${{Angiospermae_contigs[@]/$line}}" )
                    echo "$line" >>  {params.outputDir}Eukaryota_ref.fa  
                    Eukaryota_contigs=( "${{Eukaryota_contigs[@]/$line}}" )
                    switchcase="ANGIOSPERMAE"
                elif [[ ${{Eukaryota_contigs[@]}} =~ $line ]]; then
                    echo "$line" >>  {params.outputDir}Eukaryota_ref.fa
                    Eukaryota_contigs=( "${{Eukaryota_contigs[@]/$line}}" )
                    switchcase="EUKARYOTA"
                elif [[ ${{Bacteria_contigs[@]}} =~ $line ]]; then
                    echo "$line" >>  {params.outputDir}Bacteria_ref.fa  
                    Bacteria_contigs=( "${{Bacteria_contigs[@]/$line}}" )
                    switchcase="BACTERIA"
                elif [[ ${{Arbuscular_contigs[@]}} =~ $line ]]; then
                    echo "$line" >>  {params.outputDir}Arbuscular_ref.fa  
                    Arbuscular_contigs=( "${{Arbuscular_contigs[@]/$line}}" )
                    switchcase="ARBUSCULAR"
                elif [[ ${{Fungi_contigs[@]}} =~ $line ]]; then
                    echo "$line" >>  {params.outputDir}Fungi_ref.fa  
                    Fungi_contigs=( "${{Fungi_contigs[@]/$line}}" )
                    switchcase="FUNGI"
                elif [[ ${{Archaea_contigs[@]}} =~ $line ]]; then
                    echo "$line" >>  {params.outputDir}Archaea_ref.fa  
                    Archaea_contigs=( "${{Archaea_contigs[@]/$line}}" )
                    switchcase="ARCHAEA"
                elif [[ ${{Imperia_contigs[@]}} =~ $line ]]; then
                    echo "$line" >>  {params.outputDir}Imperia_ref.fa  
                    Imperia_contigs=( "${{Imperia_contigs[@]/$line}}" )
                    switchcase="IMPERIA"
                else
                    echo "$line" >>  {params.outputDir}Eukaryota_ref.fa    
                    Eukaryota_contigs=( "${{Eukaryota_contigs[@]/$line}}" )
                    switchcase="EUKARYOTA"
                fi
            else 
                case $switchcase in
                    ANGIOSPERMAE)
                        echo "$line" >>  {params.outputDir}Angiospermae_ref.fa
                        echo "$line" >>  {params.outputDir}Eukaryota_ref.fa  
                        ;;
                    EUKARYOTA)
                        echo "$line" >>  {params.outputDir}Eukaryota_ref.fa
                        ;;
                    BACTERIA) 
                        echo "$line" >>  {params.outputDir}Bacteria_ref.fa  
                        ;;
                    ARBUSCULAR)
                        echo "$line" >>  {params.outputDir}Arbuscular_ref.fa  
                        ;;
                    FUNGI)
                        echo "$line" >>  {params.outputDir}Fungi_ref.fa  
                        ;;
                    ARCHAEA)
                        echo "$line" >>  {params.outputDir}Archaea_ref.fa  
                        ;;
                    IMPERIA)
                        echo "$line" >>  {params.outputDir}Imperia_ref.fa  
                        ;;
                    *)
                        echo "$line" >>  {params.outputDir}Eukaryota_ref.fa
                        ;;
                esac
            fi
        done < {input.ref}
        """
