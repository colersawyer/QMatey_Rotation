
### change into directory holding filtered_species_sequences.txt


### This is all about done, I wanted to attempt to run the all of the qmatey data through this pipeline but I am unable to find that data on the computer that I can ssh into


### need to move files titled "Filtered_species_sequences.txt" to here.

### grab specific sequence, for this run we are using 'Metazoa' and making a new file, incorrect_reads

for i in $(ls *_filtered_species_sequences.txt);do
	grep 'Metazoa' $i > ${i%_filtered_species_sequences*}_incorrect_reads.txt 
done

###Sort incorrect reads into _sorted.txt file

for i in $(ls *_incorrect_reads.txt);do
	sort -u -t$'\t' -k2,2 $i | awk -F '\t' '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22}' > ${i%_incorrect_reads*}_sorted.txt
done

### print out hseqID as well as taxonomic information

for i in $(ls *_sorted.txt);do
	awk -F '\t' '{print $1"\t"$10"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20}' $i > ${i%_sorted*}_original_taxonomy.txt
done

### Create sequence_id.txt file by printing hseqid and ### into new file

for i in $(ls *_sorted.txt);do
	awk '{print " " $1 "\t" $8}' $i | awk 'gsub(" ", ">")' > ${i%_sorted*}_sequence_id.txt
done

### create new lines to make new_fasta files

for i in $(ls *_sequence_id.txt);do
	awk '{print $1"\n"$2}' $i > ${i%_sequence_id*}_new.fasta
done

### use fasta file for blast

export PATH=$PATH:/home/cole/Desktop/qmatey/ncbi-blast-2.9.0+/bin

export BLASTDB=/home/cole/Desktop/qmatey/ncbi-blast-2.9.0+/blast/db

for i in $(ls *_new.fasta);do ncbi-blast-2.9.0+/bin/blastn -task megablast -remote -query $i -db nt -max_target_seqs 5 -max_hsps 1 -outfmt \
"6 qseqid sseqid length mismatch evalue pident qcovs qseq sseq staxids stitle" -out ${i%_new*}_new_blast_sorted.megablast; done

### after blasting, take new blast and create new taxonomy file.

for i in $(ls *_new_blast_sorted.megablast);do
	awk -F '\t' '{print $1"\t"$10}' $i >  ${i%_new_blast_sorted*}_new_taxonomy.txt
done

### concatenate original and new_taxonomy files into one known into intermediate (itm) file

for o in $(ls *_original_taxonomy.txt);do

	n=$(basename "$o" _original_taxonomy.txt); cat $o "$n"_new_taxonomy.txt > "$n"_itm.txt

done

### Print hseqid and taxid into combined taxonomy file. I moved to a new location to make downstream stuff easier

for i in $(ls *_itm.txt);do

	awk -F '\t' '{print $1"\t"$2}' $i > /home/cole/Desktop/qmatey/ncbi-blast-2.9.0+/blast/db/new_taxdump/${i%_itm*}_combined_taxonomy.txt
done

rm *_itm.txt

cd /home/cole/Desktop/qmatey/ncbi-blast-2.9.0+/blast/db/new_taxdump

### use combined taxonomy to make file for ranked lineage grab that just has the taxids

### note: I wanted to say I do the ranked lineage here on all of the taxids (new and old) to preserve the order of everything. This was a work around I found so that all of the data doesn't get jumbled.
### This may be unnecessary and there may be a better way to do this, I just couldn't figure out how.

for i in $(ls *_combined_taxonomy.txt);do

	awk -F '\t' '{print $2}' $i > ${i%_combined_taxonomy*}_ranked_lineage_prep.txt
done


#I wrote the code below (###) to grab taxaId's that do not have ranked lineage info already. Without the ranked lineage info, the length of any column past $2, will be zero. 
#This script grabs column two if the length of column 5 is = 0, meaning it grabs all items in column two without ranked lineage info. I do not think that this is necessary, but wanted to be able to use it to save time on ranked lineage pulls.
###awk -F '\t' '!length($5) {print $2}' combined_taxonomy.txt > /home/cole/Desktop/qmatey/ncbi-blast-2.9.0+/blast/db/new_taxdump/file_for_ranked_lineage.txt

### Create tab delimited file for ranked lineage that can be use to grab info and place next to original hseqid values


awk '{gsub(/\t\t/,"\tNA\t"); print}' rankedlineage.dmp | awk '{gsub(/[|]/,""); print}' | awk '{gsub(/\t\t/,"\t"); print}' > rankedlineage_tabdelimited.dmp

echo $'tax_id\ttaxname\tspecies\tgenus\tfamily\torder\tclass\tphylum\tkingdom\tsuperkingdom\t' | cat - rankedlineage_tabdelimited.dmp > rankedlineage_edited.dmp

rm rankedlineage_tabdelimited.dmp

for i in $(ls *_ranked_lineage_prep.txt);do

	awk '{gsub(";","\n"); print}' $i | sed -e '1s/staxids/tax_id/' > ${i%_ranked_lineage_prep*}_ranked_lineage.txt

	done

rm *_ranked_lineage_prep.txt

for i in $(ls *_ranked_lineage.txt);do

	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' /home/cole/Desktop/qmatey/ncbi-blast-2.9.0+/blast/db/new_taxdump/rankedlineage_edited.dmp $i | \

	awk '{gsub(/ /,"_"); print }' > ${i%_ranked_lineage*}_unique_ranked_lineages.txt 

done

### CANT ACTUALLY REMEMBER WHAT THIS DOES NEED TO RETRY IT.

for i in $(ls *_combined_taxonomy.txt);do

	awk '!($2="")' $i > ${i%_combined_taxonomy*}_prep_for_comparison.txt

done

### Create a temporary file by concatenating the ranked lineage info with the appropriate hseqid info into  a single file

for o in $(ls *_prep_for_comparison.txt);do

	n=$(basename "$o" _prep_for_comparison.txt); paste $o "$n"_unique_ranked_lineages.txt > ${o%_prep_for_comparison*}_almost_there_file.txt

done


rm *_prep_for_comparison.txt

### Remove all of the synthetic chromosomes. This step is optional, I just wanted to get rid of them as they are a bit annoying.

for i in $(ls *_almost_there_file.txt);do

	awk '$11 != "NA"' $i > ${i%_almost_there_file*}_final_comparison_file.txt

done

### Make a new directory for all of the different datasets so that they don't get mix-matched together.

mkdir analysis

mv *_final_comparison_file.txt analysis

cd analysis

### Make Directory for individual datasets

for i in $(ls *_final_comparison_file.txt);do

	mkdir ${i%_final_comparison_file*}_joined_data 

done

### Move final comparison final into their respective folders

for i in $(ls *_final_comparison_file.txt);do

	mv $i -t ${i%_final_comparison_file*}_joined_data

done

### go into each individual folder and print out files for each hseqid. Each file will be named after the hseqid and will contain all of the new and old taxa information side by side for comparison

for d in ./*/;do

	(cd "$d" && awk '{print>$1".txt"}' *_final_comparison_file.txt && rm *_final_comparison_file.txt)

done

### go into each individual folder and remove any duplicate lines. This means that if in the second search if something only came back as itself "Homo sapiens -> Homo sapiens" then the file would
### have one line in it. This is ok though as I use any file with one line in it as a species level distinction

for d in ./*/;do

(cd "$d" && for j in $(ls *.txt);do

awk '!seen[$0]++' $j > ${j%*.txt}_sorted.txt

done)

done

##########################################################################################################################################


for d in ./*/;do

(cd "$d" && for j in $(ls *_sorted.txt);do 

awk '{print $0"\n"$0}' $j > ${j%_sorted*}_plz_delete_me.txt

done)

done

###

for d in ./*/;do

(cd "$d" && for j in $(ls *_plz_delete_me.txt);do

awk -v x=3 'NR==x{exit 1}' $j && mv $j ${j%_plz_delete_me*}_sorted.txt

done)

done

### The above code finds the right files but also deletes them, need to find away to not delete them

####COULD MOVE THEM TO A NEW FOLDER, THEN CHANGE BACK AND REMOVE OLD ONES. OR JUST REWRITE CODE TO DELETE THE CORRECT ONES

######################################################################################################################################

for d in ./*/;do

(cd "$d" && for j in $(ls *_sorted.txt);do

awk 'NR <= 1 {print; next}{print | "sort -k11,11 -r -k10,10 -r -k9,9 -r -k8,8 -r -k7,7 -r -k6,6 -r -k5,5 -r"}' $j | grep .> ${j%_sorted.*}_distance_sorting.txt

done)

done



###Take the sorted file and print a new file, where it goes hseqid, old taxa id, old tax info, new tax id, new tax info

for d in ./*/;do

(cd "$d" && for j in $(ls *_sorted.txt);do

awk 'NR == 1 {print $2"\t"$3"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11} NR == 2 {print $2"\t"$3"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' $j | tr "\n" "\t"> ${j%_sorted.*}_newboi.txt

done)

done



#Define variables suchs that 1 variable exists for each level of taxonomic information for THE ORIGINAL RANKED LINEAGE INFORMATION
#Since I sorted unique above, all of the original blasts that did not change, (i.e. they came back as themselves homo sapeins -> homo sapiens) will have only one line. Need to change the code beneath to reflect that.
#s = species
#g = genus
#f = family
#o = order
#c = class
#p = phylum
#k = kingdom
#a = super kingdom

for d in ./*/;do

	(cd "$d" &&

	for j in $(ls *_distance_sorting.txt);do

	awk 'NR == 1 {s=$3; g=$5; f=$6; o=$7; c=$8; p=$9; k=$10; a=$11}

	NR == 2 {if (s == $3) {print $1 ":\tSpecies Match"; exit} else if (g == $5) {print $1 ":\tGenus Match"; exit} else if (f == $6) {print  $1 ":\tFamily Match"; exit} \

	else if (o == $7) {print  $1 ":\tOrder Match"; exit} else if (c == $8) {print  $1 ":\tClass Match"; exit} else if (p == $9) {print  $1 ":\tPhylum Match"; exit}  \

	else if (k == $10) {print  $1 ":\tKingdom Match"; exit} else if (a == $11) {print  $1 ":\tSuper Kingdom Match"; exit} else {print "something went wrong lol"}}' $j > ${j%*_distance_sorting*}_taxonomic_shift.txt

	done)

done






#paste together the taxonomic shift with the newboi file to get all of the data in one spot. Can cat this later to get the master taxa shift 



for d in ./*/;do

(cd "$d" && 

for o in $(ls *_taxonomic_shift.txt);do

	n=$(basename "$o" _taxonomic_shift.txt); paste $o "$n"_newboi.txt > ${o%_taxonomic_shift.txt*}_taxa_shift_data.txt

done)

done


###concatenate all taxa shift files into one master taxa shift file

for d in ./*/;do

(cd "$d" && cat *_taxa_shift_data.txt > master_taxa_shift.txt)

done



###Rename all master taxa shift files according to their original dataset (tanzania dataset will generate a file called tanzania_master_taxa_shift.txt)

for d in ./*/;do

(cd "$d" && a=${PWD##*/} && mv master_taxa_shift.txt ${a}_master_taxa_shift.txt)

done

### remove all other files but the master taxa shift

for d in ./*/;do

(cd "$d" && find . ! -name "*_master_taxa_shift.txt" -type f -exec rm -f {} +)

done















