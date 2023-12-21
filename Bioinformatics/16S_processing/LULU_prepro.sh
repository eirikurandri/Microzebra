#do the pre processing neccssary for input into LULU
makeblastdb -in ASVs.fa -parse_seqids -dbtype nucl
blastn -db ASVs.fa -outfmt '6 qseqid sseqid pident' -out match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query ASVs.fa
