# Conversões do Kraken

Objetivo principal:

-   Tornar a matriz do Kraken em matriz de OTU e taxonomia, separada para:
    -   Bactérias
    -   Fungos
    -   Virus
    -   Archaea

## 1. Converter IDs em taxonomia

```shell
## Extrair a matriz BIOM
mykraken-biom \
	--max P \
	--min S \
	-o test.tsv \
	--fmt tsv \
	*.tsv
	
## Extraindo os IDs
taxonkit lineage taxids.txt | \
taxonkit reformat -f "k__{k}\tp__{p}\tc__{c}\to__{o}\tf__{f}\tg__{g}\ts__{s}" | \
cut -f1,3,4,5,6,7,8,9 > taxonomy.tsv

## Criando o taxonomy final
echo "#TAXONOMY\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies" > taxonomy_final.tsv

cat taxonomy.tsv >> taxonomy_final.tsv
```

## 2. Filtrar bactéria, fungo, vírus e archaea

Filtra da matriz de contagem!

```shell
mkdir -p filter
in="otu_table_full.txt"

cp ${in} temp_table.tsv

dirty_ids="/work/rcsilva/projects/masters/data/08_taxonIds/00_dirt2.txt"

for line in $(cat ${dirty_ids})
do
	grep -v "^${line}" temp_table.tsv > filter/f_table_tmp.txt
	mv filter/f_table_tmp.txt temp_table.tsv
done
```

### 3. Limpar IDs indesejados (planta, homo sapiens)

```shell
dirty_ids="/work/rcsilva/projects/masters/data/08_taxonIds/00_dirt2.txt"
cp otu_table_full.txt temp_table.tsv

for line in $(cat ${dirty_ids})
do
	grep -v "^${line}" temp_table.tsv > temp_f_table.tsv
	mv temp_f_table.tsv temp_table.tsv
done

mv temp_table count_table_filtered.txt
```

