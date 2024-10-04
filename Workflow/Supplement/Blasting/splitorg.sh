organism_file=Relevant_organisms.csv

Bryophyta=()
Eukaryota=()
Imperia=()
Fungi=()
Gymnospermae=() 
Pteridophyta=() 
Angiospermae=()

while IFS=, read -r Br Eu Im Fu Gy Pt An; do
    if [[ $Br = *[[:space:]]* ]]
    then
        Bryophyta+=("${Br//}")
    fi
    if [[ $Eu = *[[:space:]]* ]]
    then
        Eukaryota+=("${Eu//}")
    fi
    if [[ $Im = *[[:space:]]* ]]
    then
        Imperia+=("${Im//}")
    fi
    if [[ $Fu = *[[:space:]]* ]]
    then
        Fungi+=("${Fu//}")
    fi
    if [[ $Gy = *[[:space:]]* ]]
    then
        Gymnospermae+=("${Gy//}")
    fi
    if [[ Pt = *[[:space:]]* ]]
    then
        Pteridophyta+=("${Pt//}")
    fi
    if [[ $An = *[[:space:]]* ]]
    then
        Angiospermae+=("${An//}")
    fi
done < $organism_file


for i in ${!Imperia[@]}; do    
    echo "${Imperia[$i]}"
done