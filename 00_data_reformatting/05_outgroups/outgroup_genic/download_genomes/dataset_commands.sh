#!/bin/bash
rm -rf ncbi_dataset
rm -rf genome_sequences
rm -rf fastas 

mkdir -p genome_sequences
mkdir -p fastas

datasets download genome accession --inputfile outgroup_ids --include genome --dehydrated
unzip -o ncbi_dataset.zip
datasets rehydrate --directory .

# move all the whole genomes into a directory
find ./ncbi_dataset/ -name "*.fna" -exec mv {} ./genome_sequences/ \;

# add the anole mitochondrial genome to the whole genome assembly
cat ./EU747728.2.fasta ./genome_sequences/GCF_000090745.2_AnoCar2.0_genomic.fna > ./genome_sequences/Anole_w_mito.fa
rm ./genome_sequences/GCF_000090745.2_AnoCar2.0_genomic.fna
mv ./genome_sequences/Anole_w_mito.fa ./genome_sequences/GCF_000090745.2_AnoCar2.0_genomic.fna

# do the same with the painted turtle mitochondria
cat ./NC_002073.3.fasta ./genome_sequences/GCF_000241765.5_Chrysemys_picta_BioNano-3.0.4_genomic.fna > ./genome_sequences/Turt_w_mito.fa
rm ./genome_sequences/GCF_000241765.5_Chrysemys_picta_BioNano-3.0.4_genomic.fna
mv ./genome_sequences/Turt_w_mito.fa ./genome_sequences/GCF_000241765.5_Chrysemys_picta_BioNano-3.0.4_genomic.fna

# make genic files for each of the downloaded genomes
mito=EU747728.2
gene_tsv=../gene_tsv/Anolis_carolinensis-GCF_000090745.2.tsv
genome_file=./genome_sequences/GCF_000090745.2_AnoCar2.0_genomic.fna
pull_genic_regions.py ${gene_tsv} ${genome_file} --mito_acc ${mito}

mito=NC_002073.3
gene_tsv=../gene_tsv/Chrysemys_picta-GCF_000241765.5.tsv
genome_file=./genome_sequences/GCF_000241765.5_Chrysemys_picta_BioNano-3.0.4_genomic.fna
pull_genic_regions.py ${gene_tsv} ${genome_file} --mito_acc ${mito}

mito=NC_053523.1
gene_tsv=../gene_tsv/Gallus_gallus-GCF_016699485.2.tsv
genome_file=./genome_sequences/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna
pull_genic_regions.py ${gene_tsv} ${genome_file} --mito_acc ${mito}

mito=NC_020155.1
gene_tsv=../gene_tsv/Geotrypetes_seraphini-GCF_902459505.1.tsv
genome_file=./genome_sequences/GCF_902459505.1_aGeoSer1.1_genomic.fna
pull_genic_regions.py ${gene_tsv} ${genome_file} --mito_acc ${mito}

mito=NC_012920.1
gene_tsv=../gene_tsv/Homo_sapiens-GCF_000001405.40.tsv
genome_file=./genome_sequences/GCF_000001405.40_GRCh38.p14_genomic.fna
pull_genic_regions.py ${gene_tsv} ${genome_file} --mito_acc ${mito}

mito=NC_001804.1 
gene_tsv=../gene_tsv/Latimeria_chalumnae-GCF_000225785.1.tsv
genome_file=./genome_sequences/GCF_000225785.1_LatCha1_genomic.fna
pull_genic_regions.py ${gene_tsv} ${genome_file} --mito_acc ${mito}

mito=NC_023515.1
gene_tsv=../gene_tsv/Microcaecilia_unicolor-GCF_901765095.1.tsv
genome_file=./genome_sequences/GCF_901765095.1_aMicUni1.1_genomic.fna
pull_genic_regions.py ${gene_tsv} ${genome_file} --mito_acc ${mito}

mito=NC_006303.1
gene_tsv=../gene_tsv/Rhinatrema_bivittatum-GCF_901001135.1.tsv
genome_file=./genome_sequences/GCF_901001135.1_aRhiBiv1.1_genomic.fna
pull_genic_regions.py ${gene_tsv} ${genome_file} --mito_acc ${mito}

mito=NC_006839.1
gene_tsv=../gene_tsv/Xenopus_tropicalis-GCF_000004195.4.tsv
genome_file=./genome_sequences/GCF_000004195.4_UCB_Xtro_10.0_genomic.fna
pull_genic_regions.py ${gene_tsv} ${genome_file} --mito_acc ${mito}

mv *-genic_regions.fasta ./fastas/ 
mv -r ./fastas/ ../

# cleanup
#rm -r genome_sequences/ ncbi_dataset/ README.md
