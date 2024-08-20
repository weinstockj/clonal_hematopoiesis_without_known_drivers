in_date="2022_08_29"
out_date="2022_08_29"
threads=4
$HOME/downloads/ensembl-vep/vep --cache -i output/vcf/training_vcf_$in_date.vcf -o output/txt/training_vcf_annotated_$out_date.txt --tab --regulatory --canonical --biotype --symbol --check_existing --plugin CADD,$HOME/downloads/cadd_data/whole_genome_SNVs.tsv.gz --flag_pick --fork $threads --force_overwrite

