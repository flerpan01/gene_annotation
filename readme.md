Generates a table of the annotated genes for mouse genome (GRCm39) from the ensembl database. Using the R-packages [`biomaRt`](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) to pull gene data. Table consists of 10 columns:

1. `chr` Chromosome ID.
2. `start` Gene start.
3. `end` Gene end.
4. `strand` Strand info, given in `+`, `-` or `*` forward, reverse or no info, respectively.
5. `gene_id` Ensembl gene ID
6. `gene_name` External gene name
7. `entrez_id` Entrez ID
8. `gene_info` Gene discription
9. `gene_type` Protein coding, non-coding RNA, pseudogene etc
10. `gene_type2` A concatination of the gene types where the micro RNAs are label ncRNA and all pseudogene variants are labeled pseudogenes
