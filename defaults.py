#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Alexander Stewart"


class ProgDefaults:
    ENTREZ_EMAIL = "FILL"
    ENTREZ_API = "FILL"

    ORDER_DICT = {"species": 0, "accession": 1, "range": 2, "gene": 3, "date": 4, "length": 5}
    OUTDIR_PREFIX = "scraped_sequences"
    COMB_FILE = "all_chunks.fasta"
    CHUNK_SIZE = 2000
    EXCLUDE_TERMS = ["unisexual", "cf.", " x "]
    FOCAL_SEQFILE_NAME = "focal_scraped_sequences.fasta"
    SPP_FILE_NAME = "captured_caudate_taxa"

    MITO_GOIS = {"nadh dehydrogenase subunit 4": "ND4",
                 "nadh dehydrogenase subunit iv": "ND4",

                 "nadh dehydrogenase subunit 2": "ND2",
                 "nadh dehydrogenase subunit ii": "ND2",

                 "cytochrome c oxidase subunit 1": "CO1",
                 "cytochrome c oxidase subunit i": "CO1",
                 "cytochrome oxidase subunit 1": "CO1",
                 "cytochrome oxidase subunit i": "CO1",
                 "chytochrome c oxidase subunit 1": "CO1",
                 "chytochrome c oxidase subunit i": "CO1",

                 "cytochrome c oxidase subunit 2": "CO2",
                 "cytochrome c oxidase subunit ii": "CO2",
                 "cytochrome oxidase subunit 2": "CO2",
                 "cytochrome oxidase subunit ii": "CO2",
                 "chytochrome c oxidase subunit 2": "CO2",
                 "chytochrome c oxidase subunit ii": "CO2",

                 "12s ribosomal rna": "12S",
                 "s-rrna": "12S",
                 "small subunit ribosomal rna": "12S",
                 "12s small subunit ribosomal rna": "12S",

                 "16s ribosomal rna": "16S",
                 "16s ribosoaml rna": "16S",
                 "l-rrna": "16S",
                 "large subunit ribosomal rna": "16S",
                 "16s large subunit ribosomal rna": "16S",


                 "cytochrome b": "cytB",
                 "cytochrome b apoenzyme": "cytB",

                 "brain derived neurotrophic factor": "BDNF",
                 "brain-derived neurotrophic factor": "BDNF",

                 "rag1": "RAG1",
                 "recombinase activating protein": "RAG1",
                 "recombinase activating protein 1": "RAG1",
                 "recombination activating gene 1": "RAG1",
                 "recombination activating protein 1": "RAG1",
                 "recombination activating protein-1": "RAG1",
                 "recombination activation protein 1": "RAG1",
                 "v(d)j recombinase subunit": "RAG1",

                 "pro-opimelanocortin": "POMC",
                 "pro-opiomelanocortin": "POMC",
                 "proopiomelanocortin": "POMC",

                 "solute carrier family 8 member 3": "SLC8A3"}