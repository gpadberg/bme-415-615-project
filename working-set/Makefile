# Makefile for gene processing pipeline with copying full folders

.PHONY: all deseq2 gene_subsets get_gene_ids copy_common

all: deseq2 gene_subsets get_gene_ids copy_common

deseq2:
	python3 deseq2-to-csv.py

gene_subsets:
	python3 gene_subsets.py

get_gene_ids:
	python3 get_gene_ids.py

copy_common:
	@NUM=$$(tail -n +2 sorted_gene_subsets/both_stress_genes.csv | wc -l); \
	echo "Both-stress gene count: $$NUM"; \
	DEST="$${NUM}_common_genes"; \
	echo "Creating folder $$DEST"; \
	mkdir -p $$DEST; \
	echo "Copying folders from current directory into $$DEST"; \
	cp -r significant_gene_ids $$DEST/; \
	cp -r sorted_gene_subsets $$DEST/; \
	echo "Done: copied significant_gene_ids and sorted_gene_subsets into $$DEST"
