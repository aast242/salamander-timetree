blast_out=blastn_out6/Latimeria_chalumnae-GCF_000225785.1-genic_regions.fasta.out6
subj_fasta=../outgroup_genic/fastas/Latimeria_chalumnae-GCF_000225785.1-genic_regions.fasta
query_fasta=pyron_williams_hime_genbank.fasta
pull_gene_from_genome.py ${blast_out} ${subj_fasta} ${query_fasta} "Latimeria chalumnae"

blast_out=blastn_out6/Anolis_carolinensis-GCF_000090745.2-genic_regions.fasta.out6
subj_fasta=../outgroup_genic/fastas/Anolis_carolinensis-GCF_000090745.2-genic_regions.fasta
pull_gene_from_genome.py ${blast_out} ${subj_fasta} ${query_fasta} "Anolis carolinensis"

blast_out=blastn_out6/Chrysemys_picta-GCF_000241765.5-genic_regions.fasta.out6
subj_fasta=../outgroup_genic/fastas/Chrysemys_picta-GCF_000241765.5-genic_regions.fasta
pull_gene_from_genome.py ${blast_out} ${subj_fasta} ${query_fasta} "Chrysemys picta"

blast_out=blastn_out6/Gallus_gallus-GCF_016699485.2-genic_regions.fasta.out6
subj_fasta=../outgroup_genic/fastas/Gallus_gallus-GCF_016699485.2-genic_regions.fasta
pull_gene_from_genome.py ${blast_out} ${subj_fasta} ${query_fasta} "Gallus gallus"

blast_out=blastn_out6/Geotrypetes_seraphini-GCF_902459505.1-genic_regions.fasta.out6
subj_fasta=../outgroup_genic/fastas/Geotrypetes_seraphini-GCF_902459505.1-genic_regions.fasta
pull_gene_from_genome.py ${blast_out} ${subj_fasta} ${query_fasta} "Geotrypetes seraphini"

blast_out=blastn_out6/Homo_sapiens-GCF_000001405.40-genic_regions.fasta.out6
subj_fasta=../outgroup_genic/fastas/Homo_sapiens-GCF_000001405.40-genic_regions.fasta
pull_gene_from_genome.py ${blast_out} ${subj_fasta} ${query_fasta} "Homo sapiens"

blast_out=blastn_out6/Microcaecilia_unicolor-GCF_901765095.1-genic_regions.fasta.out6
subj_fasta=../outgroup_genic/fastas/Microcaecilia_unicolor-GCF_901765095.1-genic_regions.fasta
pull_gene_from_genome.py ${blast_out} ${subj_fasta} ${query_fasta} "Microcaecilia unicolor"

blast_out=blastn_out6/Rhinatrema_bivittatum-GCF_901001135.1-genic_regions.fasta.out6
subj_fasta=../outgroup_genic/fastas/Rhinatrema_bivittatum-GCF_901001135.1-genic_regions.fasta
pull_gene_from_genome.py ${blast_out} ${subj_fasta} ${query_fasta} "Rhinatrema bivittatum"

blast_out=blastn_out6/Xenopus_tropicalis-GCF_000004195.4-genic_regions.fasta.out6
subj_fasta=../outgroup_genic/fastas/Xenopus_tropicalis-GCF_000004195.4-genic_regions.fasta
pull_gene_from_genome.py ${blast_out} ${subj_fasta} ${query_fasta} "Xenopus tropicalis"
