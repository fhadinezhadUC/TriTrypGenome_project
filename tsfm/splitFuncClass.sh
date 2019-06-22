#!/bin/bash

# functions we have from homo
# A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  X  Y X
# we need to go through all the files in this path and split the functional class

#filenames=$(ls /home/fatemeh/Leishmania_2019/Leishmania_2019/Results/tsfmInput-output/input2/*.fasta)
#pathstring="/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/TriTrypgenemodel_Intersection/"
pathstring="/home/fatemeh/TriTrypGenome_project/tsfm/input2/"

filenames=$(ls "${pathstring}"*.fasta)
array=(A R N D C Q E G H I L K M X F P S T W Y V X)

for ref in $filenames; do
    ref2=${ref#"${pathstring}"}
    ref3=${ref2%".fasta"}
    cd  "${pathstring}"
    mkdir "$ref3"
    printf " %s\n" $ref3
mkdir "${pathstring}$ref3/clustalW"
for item in ${array[*]}
do
    printf "   %s\n" $item
    cat "$ref" | fasgrep "_${item}$" > "${pathstring}$ref3/${item}.fasta"
    cat "${pathstring}$ref3/${item}.fasta" | fasconvert -o clustalw > "${pathstring}$ref3/clustalW/${ref3}_${item}.aln"
  

done
done
