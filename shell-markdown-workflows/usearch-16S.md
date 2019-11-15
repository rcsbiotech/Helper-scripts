# Microbiomas com USEARCH11

Tutorial para usar o USEARCH gratuito (v.11), útil para amostras de 16S ou ITS de tamanho de amplicon **conhecido**, com dados de no máximo 3gb. 

## Alertas

-   O passo de orientação das reads `ITS` está cagando na amostragem, evitando usar por hora, até fixar as referências.

[TOC]

## Programas e dados necessários

Na `pipeline` recomendo dois programas:

- USEARCH11 (32bit para até 4.0gb de sequências unidas, 64bit para acima);
- FastQC.

O USEARCH11 de 32 bits é gratuito, e deve ser solicitado [neste link](https://drive5.com/usearch/download.html), colocando seu e-mail para cadastro único.

O FastQC também é gratuito, e pode ser baixado [neste link](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip).

Na versão deste tutorial, servem apenas dados pareados de sequenciamento Illumina.

## Preâmbulo: estrutura de pasta, informações necessárias para a análise e nome de variáveis

Para começar a fazer as análises, são necessárias as seguintes informações, que são alimentadas no programa

- Os _primers_ utilizados;
- O tamanho de _amplicon_ esperado.

São passadas as variáveis com os seguintes nomes nos processos (copie e cole isto no seu terminal cada vez que for rodar o estudo do começo ao fim):

```shell
## Número de processadores
p_threads=35

## Variáveis de nome, refaça conforme os testes que desejar rodar
merge_no=M01
prepr_no=P01
annot_no=A01

## Parâmetros do programa: escolher bacteria ou fungi
p_organism="bacteria"
annotation_parameter="bacteria"

## Parâmetros da etapa merge
p_diffs=30
p_id=90

## Parâmetro de erros esperados
p_ee=1.00

## Parâmetros de truncagem e trimagem
p_trunc=210
p_trimR1=20
p_trimR2=21

## Parâmetros de anotação
	# Escolher entre 0.0 e 1.0
	# Bactérias: 0.85 - 0.99
	# Fungos: 0.5+
p_sintax_cutoff=0.5

# Default de bactérias
annotation_database="rdp.mito.udb"
# Default de fungos
# annotation_database="rdp-its.udb"	  
```

Essas variáveis são utilizadas na nomenclatura das pastas, que irão salvar os resultados das análises de forma seriada, para não perder nenhum dado ou salvar em cima de análises já feitas.

Não bastando, eu sugiro que você sempre faça uma tabela como a abaixo para controlar os parâmetros passados para cada rodada, como da forma abaixo:

#### Tabela: união das sequências

| Merge | Diferenças | Identidade |
| ----- | ---------- | ---------- |
| M01 | 20 | 80% |

#### Tabela: pré-processamento

| Prepr_no | EE | Trunc | Trim_R1 | Trim_R2 |
| -------- | -- | ----- | ------- | ------- |
| P01 | 01 | 400 | 19 | 20 |

#### Tabela: Anotação

| Annot_no | DB | OTUS/zOTUS |
| -------- | -- | ---------- |
| A01 | RDP_small | zOTUS |

## Criação dos primeiros diretórios e movimentação dos dados

Crie uma pasta vazia, e nessa, passe as seguintes pastas e comandos:

```shell
mkdir -p ./dados/brutos
mkdir -p ./dados/oligos
mkdir -p ./dados/referencias
mkdir -p ./dados/intel
mkdir -p ./analises
mkdir -p ./analises/E01_A_fastx_info/logs
mkdir -p ./analises/E01_B_fastq_eestats2/logs
mkdir -p ./analises/E01_C_search_oligodb/logs
mkdir -p ./analises/E01_D_fastqc
mkdir -p ./analises/E02_mergepairs/logs
mkdir -p ./analises/E03_trimming/logs
mkdir -p ./analises/E04_orient/logs
mkdir -p ./analises/E05_filter/logs
mkdir -p ./analises/E06_derrep/logs
mkdir -p ./analises/E07_singletons/logs
mkdir -p ./analises/E08_pick_otus_zotus/logs
mkdir -p ./analises/E09_otu_tables/logs
mkdir -p ./analises/E10_sintax_taxonomy
mkdir -p ./analises/E11_distance_tree
mkdir -p ./analises/E12_purge_chlr_mito
mkdir -p ./analises/E13_microbiome-analyst
mkdir -p ./analises/E14_work_in_R
mkdir -p ./analises/E15_PICRUSt
mkdir -p ./analises/E16_alpha
mkdir -p ./analises/E17_beta
mkdir -p ./analises/tmp

## Servidor Thor!
ln -s /data/db/microbiome/*udb ./dados/referencias/
ln -s /data/db/microbiome/oligos/* ./dados/oligos/
ln -s /data/db/microbiome/fasta/gg97.fa ./dados/referencias/

## Servidor geral
# microbiome_dir="/work/db/microbiome"
# ln -s ${microbiome_dir}/*udb ./dados/referencias/
echo 'Okay!'

```

A seguir, os dados brutos (formato `fastq`) devem ser movidos para a pasta `./dados/brutos`, e lá eles permanecerão sem receber alterações até o final do processo. Isso garante que eles fiquem integros, e que você não tenha que baixar nem descomprimí-los novamente.

:::danger
É muito importante que você rode **TODOS** os comandos de dentro da pasta criada, isto é, aquela que contém a pasta dados e análises!
:::

## Etapa 01: Compreendendo seu sequenciamento

Com o USEARCH, Edgar insiste em uma política de olhar com cuidado no sequenciamento sem correr e apertar o botão da `pipeline`, e assim, ele ensina o que deve ser investigado no sequenciamento, buscando regiões de baixa qualidade.

Nesta primeira etapa, a ideia é entender quão bom foi seu sequenciamento, quanto dele irá se aproveitar para obter as comunidades de microrganismos, e o que deverá ser feito para isso. Então, vamos a eles:

Assim, aqui são utilizados quatro comandos:

- `fastx_info`
- `fastq_eestats2`
- `search_oligodb`
- `fastqc`

### A. `fastx_info`: Sumário do sequenciamento

Para este comando, é necessário dar como entrada os arquivos **fastq**, na mesma pasta os arquivos R1 e R2 _forward_ e _reverse_, e serão caputrados no laço _for_:

O comando é executado da seguinte forma:

```shell
## Relatório de erros esperados
## Para cada 25000 sequências

## Para cada arquivo fastq, faça
for fqfile in ./dados/brutos/*R1*fastq*
do

    ## Capture o 'nome base'
    bn=$(basename $fqfile .fastq)
    fq2=$(echo $fqfile | sed 's/R1/R2/')
    
    head -n 100000 ${fqfile} > ./analises/tmp/tmp_R1.fq
    head -n 100000 ${fq2} > ./analises/tmp/tmp_R2.fq
    
    
    ## Rode o comando para cada sequência
    usearch11 \
        -fastx_info ./analises/tmp/tmp_R1.fq \
        -output ./analises/E01_A_fastx_info/${bn}.txt \
        1> ./analises/E01_A_fastx_info/logs/${bn}.stdout \
        2> ./analises/E01_A_fastx_info/logs/${bn}.stderror
        
done
```

Feito isto, na pasta de análises, em `E01_A_fastx_info`, cada arquivo de leituras  possuirá um sumário com a seguinte formatação:

```
File size 239M, 356.1k seqs, 106.0M letters and quals
Lengths min 35, lo_quartile 301, median 301, hi_quartile 301, max 301
Letter freqs G 34.0%, A 23.8%, C 22.7%, T 19.4%, N 0.1%
0% masked (lower-case)
ASCII_BASE=33
EE mean 4.2; min 0.0, lo_quartile 2.0, median 3.3, hi_quartile 5.3, max 22.1
```

Isto exibe a média de tamanho das suas sequências por quartil, além da frequência das bases (pegar erros e conteúdo não esperado pra 16S), além da média de erros, `EE mean`.

Isto é interessante para descobrir se deu tudo certo com o sequenciamento, se as sequencias tem o tamanho correto, se o número médio de erros é compatível (menor de 05 para R1 e menor de 10 para R2).

A seguir, o próximo comando irá exibir as opções de corte por valor de erro, isto é, quanto das suas sequências você irá perder para cada valor de erro e truncagem escolhidos.

### B. `fastq_eestats2`: Estatísticas de erros do sequenciamento

Neste próximo comando, bem mais importante que o anterior, são também utilizados os arquivos fastq um por vez, e exibe o perfil de erros e cumprimento das sequências.

Idealmente, o erro máximo deve ser de $Q < 15$, o que dá 1 erro a cada 80 bases.

```shell
## Para cada arquivo fastq, faça
for fqfile in ./dados/brutos/*fastq
do

    ## Capture o 'nome base'
    bn=$(basename $fqfile .fastq)
    fq2=$(echo $fqfile | sed 's/R1/R2/')
    
    head -n 100000 ${fqfile} > ./analises/E01_B_fastq_eestats2/tmp_R1.fq
    head -n 100000 ${fq2} > ./analises/E01_B_fastq_eestats2/tmp_R2.fq
    
    ## Rode o comando para cada sequência
    usearch11 \
        -fastq_eestats2 ./analises/tmp/tmp_R1.fq \
        -output ./analises/E01_B_fastq_eestats2/${bn}.txt \
        -threads ${p_threads} \
        -length_cutoffs 150,400,20 \
        -ee_cutoffs 0.5,1.0,2.0,3.0,4.0 \
        1> ./analises/E01_B_fastq_eestats2/logs/log.${bn}.stdout \
        2> ./analises/E01_B_fastq_eestats2/logs/log.${bn}.stderr
        
done
```

### C. `search_oligodb`: Busca dos _primers_

1. Subamostrar as sequencias
2. Remover os oligos

- Objetivo: buscar os oligos nas sequencias, para posterior remoção
- `dados/oligos/oligos.fa` arquivo de oligos

```shell
for fqfile in ./dados/brutos/*R1*fastq
do

    ## Capture o 'nome base'
    bn=$(basename $fqfile .fastq)
    
    usearch11 \
        -search_oligodb $fqfile \
        -db /data/db/microbiome/oligos/oligos.fa \
        -strand both \
        -userout ./analises/E01_C_search_oligodb/${bn}.report.txt \
        -userfields query+target+qstrand+diffs+tlo+thi+trowdots \
        1> ./analises/E01_C_search_oligodb/logs/log.${bn}.stdout \
        2> ./analises/E01_C_search_oligodb/logs/log.${bn}.stderr

done
```

### D. `fastqc`: Análise visual das sequências

```shell
## Para cada arquivo fastq, faça
for fqfile in ./dados/brutos/*fastq
do

	## Cria o diretório para cada sequência
	mkdir -p ./analises/E01_D_fastqc/${fqfile}
    
    ## Gere o relatório visual
    fastqc ${fqfile} \
    	-o ./analises/E01_D_fastqc/${fqfile} \
    	--adapters ./dados/oligos/oligos.fa
        
done
```



## Etapa 02: Unindo os pares, `fastq_mergepairs`


```shell
## Default
    usearch11 \
        -fastq_mergepairs ./dados/brutos/*R1*.fastq \
        -fastqout analises/E02_mergepairs/merged.${merge_no}.fq \
        -relabel @ \
        -fastq_maxdiffs ${p_diffs} \
        -fastq_pctid ${p_id} \
        -report analises/E02_mergepairs/logs/merge.${merge_no}.report.txt \
            1> ./analises/E02_mergepairs/logs/${merge_no}.stdout \
            2> ./analises/E02_mergepairs/logs/${merge_no}.stderror
	
## Usearch sem logs
usearch11 \
	-fastq_mergepairs ./dados/brutos/*R1*.fastq \
	-fastqout analises/E02_mergepairs/merged.${merge_no}.fq \
	-relabel @ \
	-fastq_maxdiffs ${p_diffs} \
	-fastq_pctid ${p_id} \
	-report analises/E02_mergepairs/logs/merge.${merge_no}.report.txt
        
## USEARCH 8
usearch81 \
	-fastq_mergepairs dados/brutos/*R1*.fastq \
	-fastqout analises/E02_mergepairs/merged.${merge_no}.fq \
	-relabel @ \
    1> ./analises/E02_mergepairs/logs/${merge_no}.stdout \
    2> ./analises/E02_mergepairs/logs/${merge_no}.stderror
    
usearch81 \
	-fastq_mergepairs dados/brutos/*R1*.fastq \
	-fastqout analises/E02_mergepairs/merged.${merge_no}.fq \
	-relabel @

## Sequências com cabeçalhos estourados (i.e. bibliotecas perdidas)
```

### Experimental - PEAR para seqs .fudidas

```shell
# Cria o diretório
mkdir -p ./analises/E02_mergepairs/pear
mkdir -p ./analises/E02_mergepairs/relabeled
mkdir -p ./analises/E02_mergepairs/clean

p_overlap=05

## Roda o PEAR para cada um
for fq1 in $(ls ./dados/brutos/*R1*.fastq)
do
	
	## Se o arquivo não existe, então...
	
	if [[ ! -f ./analises/E02_mergepairs/pear/${fqb} ]]; then
		fq2=$(echo ${fq1} | sed 's/R1/R2/')
		fqb=$(basename ${fq1} .fastq)
	
		pear \
			--forward-fastq ${fq1} \
			--reverse-fastq ${fq2} \
			--threads 50 \
			--min-overlap ${p_overlap} \
			--output ./analises/E02_mergepairs/pear/${fqb}
		
	fi
		
	fqb_clean=$(basename ${fqb} .assembled.fastq)
	sample_id=$(echo ${fqb_clean} | sed 's/\(G[0-9]\+\).*/\1/')
				
	echo
	echo
	echo "ID: ${sample_id}"
	echo "Basename: ${fqb_clean}"
	sleep 0
	echo
	echo
				
	usearch11 \
    	-fastq_filter \
    		./analises/E02_mergepairs/pear/${fqb}.assembled.fastq \
    	-fastqout \
    		./analises/E02_mergepairs/clean/${fqb_clean}.clean.fq \
    	-relabel ${sample_id}.
    		
done
```



## Etapa 03: Truncagem e trimagem, `fastx_truncate` 

```SHELL
## Até 4.0gb
usearch11 \
	-fastx_truncate analises/E02_mergepairs/merged.${merge_no}.fq \
	-trunclen ${p_trunc} \
	-stripleft ${p_trimR1} \
	-stripright ${p_trimR2} \
	-fastqout ./analises/E03_trimming/trimmed.${merge_no}.${p_trunc}.fq
	
## Acima de 4.0gb - só trunca
usearch81 \
	-fastx_truncate analises/E02_mergepairs/merged.${merge_no}.fq \
	-trunclen ${p_trunc} \
	-fastqout ./analises/E03_trimming/trimmed.${merge_no}.${p_trunc}.tmp.fq

## Acima de 4.0gb - tira as pontas
usearch81 \
	-fastx_truncate ./analises/E03_trimming/trimmed.${merge_no}.${p_trunc}.tmp.fq \
	-stripleft ${p_trimR1} \
	-stripright ${p_trimR2} \
	-fastqout ./analises/E03_trimming/trimmed.${merge_no}.${p_trunc}.fq
	
## Acima de 4gb single end
usearch81 \
	-fastx_truncate ./analises/E03_trimming/trimmed.${merge_no}.${p_trunc}.tmp.fq \
	-stripleft ${p_trimR1} \
	-fastqout ./analises/E03_trimming/trimmed.${merge_no}.${p_trunc}.fq
	
echo "Trimming and truncating OK"
```

## Etapa 04: Orientando as leituras, `orient`

Nem ando fazendo, por não observar nenhuma vantagem prática.

```SHELL
## Orienta em relação ao banco de dados
usearch11 \
	-orient ./analises/E03_trimming/trimmed.${merge_no}.${p_trunc}.fq \
	-db ./dados/referencias/gg97.udb \
	-fastqout ./analises/E04_orient/oriented.${merge_no}.${p_trunc}.fq \
	-tabbedout ./analises/E04_orient/orientation.${merge_no}.${p_trunc}.tsv
	
## Para fungos
usearch11 \
	-orient ./analises/E03_trimming/trimmed.${merge_no}.${p_trunc}.fq \
	-db ./dados/referencias/rdp-its.udb \
	-fastqout ./analises/E04_orient/oriented.${merge_no}.${p_trunc}.fq \
	-tabbedout ./analises/E04_orient/orientation.${merge_no}.${p_trunc}.tsv
	
## Indisponível para acima de 4gb!
```

## Etapa 05: Filtragem por qualidade, `fastq_filter`

```shell
## Menos que 4.0 gb
usearch11 \
	-fastq_filter ./analises/E04_orient/oriented.${merge_no}.${p_trunc}.fq \
	-fastq_maxee ${p_ee} \
	-fastaout ./analises/E05_filter/filtered.${prepr_no}.EE_${p_ee}.fa
	
## Mais que 4.0 gb
usearch81 \
	-fastq_filter ./analises/E03_trimming/trimmed.${merge_no}.${p_trunc}.fq \
	-fastq_maxee ${p_ee} \
	-fastaout ./analises/E05_filter/filtered.${prepr_no}.EE_${p_ee}.fa
```

## Etapa 06: Desreplicação, `fastx_uniques`

```shell
# 32bit (v11)
usearch11 \
	-fastx_uniques ./analises/E05_filter/filtered.${prepr_no}.EE_${p_ee}.fa \
	-fastaout ./analises/E06_derrep/derrep.${prepr_no}.EE_${p_ee}.fa \
	-sizeout \
	-relabel Uniq
    
# 64bit (v81)
usearch81 \
	-derep_fulllength ./analises/E05_filter/filtered.${prepr_no}.EE_${p_ee}.fa \
	-fastaout ./analises/E06_derrep/derrep.${prepr_no}.EE_${p_ee}.fa \
	-sizeout \
	-relabel Uniq
```

## Etapa 07: Remoção de sequências únicas, `sortbysize`

```shell
# Mínimo de 02 sequências
usearch11 \
	-sortbysize ./analises/E06_derrep/derrep.${prepr_no}.EE_${p_ee}.fa \
	-fastaout ./analises/E07_singletons/nosingle.${prepr_no}.EE_${p_ee}.fa \
	-minsize 2
```

## Etapa 08: Montagem das OTUs, `cluster_otus` e `unoise3`

#### `cluster_otus`, 97% de identidade

```shell
## Cluster a 97% - considerado obsoleto
# Também detecta quimeras com UCHIME.
usearch11 \
	-cluster_otus ./analises/E07_singletons/nosingle.${prepr_no}.EE_${p_ee}.fa \
	-otus ./analises/E08_pick_otus_zotus/otus.${prepr_no}.fa \
	-relabel Otu 
```

#### `unoise3`

```shell
## zOTUs
usearch11 \
	-unoise3 ./analises/E07_singletons/nosingle.${prepr_no}.EE_${p_ee}.fa \
	-zotus ./analises/E08_pick_otus_zotus/zotus.${prepr_no}.fa
```

## Etapa 09: Construção da tabela de OTUs, `otutab`

#### Sequências trimmadas sem reorientação

```shell
## Para OTUs clusterizadas sem passo de reorientação das leituras
usearch11 \
	-otutab ./analises/E03_trimming/trimmed.${merge_no}.${p_trunc}.fq \
	-id 0.99 \
	-otus ./analises/E08_pick_otus_zotus/otus.${prepr_no}.fa \
	-otutabout ./analises/E09_otu_tables/otutab.${prepr_no}.tsv \
	-mapout ./analises/E09_otu_tables/mapout.${prepr_no}.tsv
```

#### Outros

```shell
## Para OTUs clusterizadas
usearch11 \
	-otutab ./analises/E04_orient/oriented.${merge_no}.${p_trunc}.fq \
	-id 1.00 \
	-otus ./analises/E08_pick_otus_zotus/otus.${prepr_no}.fa \
	-otutabout ./analises/E09_otu_tables/otutab.${prepr_no}.tsv \
	-mapout ./analises/E09_otu_tables/mapout.${prepr_no}.tsv
	

	
## Para zOTUs
usearch11 \
	-otutab ./analises/E03_trimming/trimmed.${merge_no}.${p_trunc}.fq \
	-otus ./analises/E08_pick_otus_zotus/zotus.${prepr_no}.fa \
	-otutabout ./analises/E09_otu_tables/zotutab.${prepr_no}.tsv \
	-id 1.00 \
	-mapout ./analises/E09_otu_tables/zmapout.${prepr_no}.tsv
	
# (OK): Caso não tenha sido orientado (64bits)
# (OK): vsearch (64bit)
fq2fa ./analises/E03_trimming/trimmed.${merge_no}.${p_trunc}.fq ./analises/E03_trimming/trimmed.${merge_no}.${p_trunc}.fa

vsearch --usearch_global \
	./analises/E03_trimming/trimmed.${merge_no}.${p_trunc}.fa \
	--db ./analises/E08_pick_otus_zotus/otus.${prepr_no}.fa \
	--otutabout ./analises/E09_otu_tables/otutab.${prepr_no}.tsv \
	--id 1.00
	
echo "OTUS ok!"


## VSEARCH_GLOBAL para zOTUs
vsearch --usearch_global \
	./analises/E03_trimming/trimmed.${merge_no}.${p_trunc}.fa \
	--db ./analises/E08_pick_otus_zotus/zotus.${prepr_no}.fa \
	--otutabout ./analises/E09_otu_tables/zotutab.${prepr_no}.tsv \
	--id 1.00
	
echo "OTUS ok!"
```

## Etapa 10: Taxonomia `sintax`

```shell
# OTUs
# RDP database large: rdp-large.udb (bacteria)
# silva isolates: silva_ltp.udb (bacteria)
# RDP mais mitocôndrias: rdp.mito.udb (bacteria)
# Greengenes (não-recomendado): ggLarge.udb | gglarge.vsearch.udb
# Silva (não-recomendado): silva-all-seqs.udb

## Bancos de dados de fungos ##
# rdp-its.udb: fungi
# utax_2019.udb: fungi

# Escolher novamente o banco?
annotation_database="rdp.mito.udb"
p_sintax_cutoff=0.9

usearch11 \
	-sintax ./analises/E08_pick_otus_zotus/otus.${prepr_no}.fa \
	-db /data/db/microbiome/${annotation_database} \
	-tabbedout ./analises/E10_sintax_taxonomy/annot.${prepr_no}.${annot_no}.txt \
	-strand both \
	-sintax_cutoff ${p_sintax_cutoff}
	
## Remoção de sujeiras da anotação
## Trocar cabeçalhos K por P
sed -i 's/k:Bacteria/d:Bacteria/g' ./analises/E10_sintax_taxonomy/annot.${prepr_no}.${annot_no}.txt

sed -i 's/"//g' ./analises/E10_sintax_taxonomy/annot.${prepr_no}.${annot_no}.txt


```

#### Deconvolução de *taxa* (deparse_taxa.R)

```shell
# Converte a tabela do formato canônico, para o formato sem OTUs repetitivas
deparse_taxa.R \
	--input=analises/E10_sintax_taxonomy/annot.${prepr_no}.${annot_no}.txt \
	--out=analises/E10_sintax_taxonomy/annot.${annot_no}.deparsed.txt
```

#### Para zOTUs

```shell
annot_no="A02"

# ZOTUs
usearch11 \
	-sintax ./analises/E08_pick_otus_zotus/zotus.${prepr_no}.fa \
	-db ./dados/referencias/rdp.mito.udb \
	-tabbedout ./analises/E10_sintax_taxonomy/zannot.${prepr_no}.${annot_no}.txt \
	-strand both \
	-sintax_cutoff ${p_sintax_cutoff}
    
sed -i 's/k:Bacteria/p:Bacteria/g' ./analises/E10_sintax_taxonomy/annot.${prepr_no}.${annot_no}.txt

sed -i 's/k:Archaea/p:Archaea/g' ./analises/E10_sintax_taxonomy/annot.${prepr_no}.${annot_no}.txt
```

### Etapa 10.1: `SINTAX` do `VSEARCH` para acima de 4gb (64bit)

Única alternativa para rodar dados maiores que 4gb // bancos maiores que 4gb

```shell
annot_no="T01"

## Otus
vsearch10 \
	--sintax ./analises/E08_pick_otus_zotus/zotus.${prepr_no}.fa \
	--db /work/db/microbiome/gglarge.vsearch.udb \
	--tabbedout ./analises/E10_sintax_taxonomy/annot.${prepr_no}.${annot_no}.txt \
	--strand both \
	--sintax_cutoff 0.9

vsearch10 \
	--sintax ./analises/E08_pick_otus_zotus/zotus.${prepr_no}.fa \
	--db /work/db/microbiome/silva.132.vsearch.udb \
	--tabbedout ./analises/E17_Silva/annot.${prepr_no}.${annot_no}.txt \
	--strand both \
	--sintax_cutoff 0.9
	
## Para fazer um banco SINTAX
```

## Etapa 11.1: Construção de árvore de distância

#### Árvore para tabela de OTUs

```shell
## I. OTUs
usearch11 \
	-calc_distmx ./analises/E08_pick_otus_zotus/otus.${prepr_no}.fa \
	-tabbedout ./analises/E11_distance_tree/${prepr_no}_otus_mx.txt \
	-maxdist 0.2 \
	-termdist 0.3

## 2. Cria a árvore newick
usearch11 \
	-cluster_aggd ./analises/E11_distance_tree/${prepr_no}_otus_mx.txt \
    -treeout ./analises/E11_distance_tree/${prepr_no}.otus.nwk \
    -clusterout ./analises/E11_distance_tree/${prepr_no}.clusters.txt \
    -id 0.80 \
    -linkage min
```

#### Outros

```shell
# 1. Calcular a matriz de Distância

## II. zOTUs
usearch11 \
	-calc_distmx ./analises/E08_pick_otus_zotus/zotus.${prepr_no}.fa \
	-tabbedout ./analises/E11_distance_tree/${prepr_no}_zotus_mx.txt \
	-maxdist 0.2 \
	-termdist 0.3
	   
## zOTUs
usearch11 \
	-cluster_aggd ./analises/E11_distance_tree/${prepr_no}_zotus_mx.txt \
    -treeout ./analises/E11_distance_tree/${prepr_no}.zotus.nwk \
    -clusterout ./analises/E11_distance_tree/${prepr_no}.zclusters.txt \
    -id 0.80 \
    -linkage min
    
    
## USEARCH8 ##
usearch81 \
	-cluster_aggd ./analises/E11_distance_tree/${prepr_no}_otus_mx.txt \
    -treeout ./analises/E11_distance_tree/${prepr_no}.otus.nwk \
    -clusterout ./analises/E11_distance_tree/${prepr_no}.clusters.txt \
    -id 0.80 \
    -linkage min
```

## Etapa 11.2: Análises de diversidade no MicrobiomeAnalyst

```shell
## 3 arquivos necessários:
# 1. OTU table
# 2. Taxonomia
# 3. Árvore newick
annotation_parameter=bacteria

## 1b. OTU table
sed 's/OTU ID/NAME/g' ./analises/E09_otu_tables/zotutab.${prepr_no}.tsv > \
	./analises/E13_microbiome-analyst/${prepr_no}_zotu_table.txt
	
## 2. Taxonomia
parseSintax \
	./analises/E10_sintax_taxonomy/annot.${prepr_no}.${annot_no}.txt \
	${annotation_parameter} \
	./analises/E13_microbiome-analyst/${prepr_no}.${annot_no}_taxonomy.txt

## 3. Newick tree
cp ./analises/E11_distance_tree/${prepr_no}.otus.nwk ./analises/E13_microbiome-analyst/


```

#### Para zOTUs

```shell
# [ .. Tabelas de zOTUs .. ] #
	
## 1b. zOTU table
sed 's/OTU ID/NAME/g' ./analises/E09_otu_tables/zotutab.${prepr_no}.tsv > \
	./analises/E13_microbiome-analyst/${prepr_no}_zotu_table.txt

## 2. Taxonomia
parseSintax \
	./analises/E10_sintax_taxonomy/zannot.${prepr_no}.${annot_no}.txt \
	fungi \
	./analises/E13_microbiome-analyst/${prepr_no}_ztaxonomy.txt
	
## 3. Newick tree
cp ./analises/E11_distance_tree/${prepr_no}.zotus.nwk ./analises/E13_microbiome-analyst/
```

### Anotação com Greengenes

```shell
# OTUs
# RDP database large: rdp-large.udb (bacteria)
# silva isolates: silva_ltp.udb (bacteria)
# RDP mais mitocôndrias: rdp.mito.udb (bacteria)
# Greengenes (não-recomendado): ggLarge.udb | gglarge.vsearch.udb
# Silva (não-recomendado): silva-all-seqs.udb

## Bancos de dados de fungos ##
# rdp-its.udb: fungi
# utax_2019.udb: fungi

# Escolher novamente o banco?
annotation_database="gglarge.vsearch.udb"
annot_no=A03
p_sintax_cutoff=0.9

vsearch10 \
	--sintax ./analises/E08_pick_otus_zotus/zotus.${prepr_no}.fa \
	--db /work/db/microbiome/${annotation_database} \
	--tabbedout ./analises/E10_sintax_taxonomy/gg.${prepr_no}.${annot_no}.txt \
	--strand both \
	--sintax_cutoff 0.9
	
## Remoção de sujeiras da anotação
## Trocar cabeçalhos K por P
sed -i 's/p:Bacteria/d:Bacteria/g' ./analises/E10_sintax_taxonomy/annot.${prepr_no}.${annot_no}.txt

sed -i 's/"//g' ./analises/E10_sintax_taxonomy/annot.${prepr_no}.${annot_no}.txt
```



## E12. Remover cloroplastos/mitocôndrias do MicrobiomeAnalyst

```shell
## Cria a pasta de limpeza
mkdir -p ./analises/E13_microbiome-analyst/cleanse
mkdir -p ./analises/E13_microbiome-analyst/parse
cleanse_dir="./analises/E13_microbiome-analyst/cleanse"

### Copia a tabela de OTUs
cp ./analises/E13_microbiome-analyst/${prepr_no}_otu_table.txt ${cleanse_dir}

### Captura as linhas sujas
grep -E "[Cc]hloroplast" ./analises/E13_microbiome-analyst/${prepr_no}.${annot_no}_taxonomy.txt > ${cleanse_dir}/lines.txt

grep -E "[Mm]itochondria" ./analises/E13_microbiome-analyst/${prepr_no}.${annot_no}_taxonomy.txt >> ${cleanse_dir}/lines.txt

## Separa as Otus
cut -f1 ${cleanse_dir}/lines.txt > ${cleanse_dir}/headers.txt

## Separa todas as linhas que não são bagunça
for otu in $(cat ${cleanse_dir}/headers.txt); do
	echo ${otu}
	grep -Ev "${otu}\s" ${cleanse_dir}/${prepr_no}_otu_table.txt >> ${cleanse_dir}/tmp
	mv ${cleanse_dir}/tmp ${cleanse_dir}/${prepr_no}_otu_table.txt
done

## Copia o arquivo de volta para a pastas do MA
cp ${cleanse_dir}/${prepr_no}_otu_table.txt ./analises/E13_microbiome-analyst/${prepr_no}_otu_table_clean.txt
```

## Adendo I: Conversão da tabela BIOM para organismos picados

```shell
# datadir deve conter:
# otu_table.txt - matriz de OTUs
# sintax.microbiomeAnalyst.txt - arquivo de taxonomia formatado
# metadata.tsv - arquivo de metadados do QIIME2
source activate qiime2

## Copia a taxonomia
cp ./analises/E13_microbiome-analyst/${prepr_no}.${annot_no}_taxonomy.txt ./analises/E14_work_in_R/sintax.microbiomeAnalyst.txt

## Copia a OTU-table
cp ./analises/E13_microbiome-analyst/${prepr_no}_otu_table.txt ./analises/E14_work_in_R/otu_table.txt

## Copia os metadados
cp ./dados/intel/metadados.tsv ./analises/E14_work_in_R/metadata.txt

datadir="./analises/E14_work_in_R"
processed="${datadir}/processed"
output="${datadir}/output"
mkdir -p ${output}
mkdir -p ./tmp
 
## Teste: parse
sed 's/\t/;/g' ${datadir}/sintax.microbiomeAnalyst.txt | \
sed 's/^\(Otu[0-9]*\);/\1\t/g' | \
sed 's/$/\t1.00/g' |
sed 1d > ${datadir}/tax.tmp
sed -i '1 i\#OTUID\ttaxonomy\tconfidence' ${datadir}/tax.tmp

## Remoção de cloroplastos
grep -E "[cC]hloroplast" ${datadir}/tax.tmp > ${datadir}/chloro.lines.txt
grep -oh "Otu[0-9]*" ${datadir}/chloro.lines.txt > ${datadir}/chloro.columns.txt

for chline in $(cat ${datadir}/chloro.columns.txt)
do
	grep -v "${chline}" ${datadir}/otu_table.txt > ${datadir}/otu_tmp.txt
	mv ${datadir}/otu_tmp.txt ${datadir}/otu_table.txt
done

## Remoção de mitocôndrias

grep -E "[mM]itochondria" ${datadir}/tax.tmp > ${datadir}/mito.lines.txt
grep -oh "Otu[0-9]*" ${datadir}/mito.lines.txt > ${datadir}/mito.columns.txt

for mtline in $(cat ${datadir}/mito.columns.txt)
do
	grep -v "${mtline}" ${datadir}/otu_table.txt > ${datadir}/otu_tmp.txt
	mv ${datadir}/otu_tmp.txt ${datadir}/otu_table.txt
done

## Remoção de ruídos numéricos
grep -E "1\s1\s1\s1" ${datadir}/tax.tmp > ${datadir}/number.lines.txt
grep -E "2\s2\s2\s2" ${datadir}/tax.tmp > ${datadir}/number.lines.txt
grep -E "3\s3\s3\s3" ${datadir}/tax.tmp > ${datadir}/number.lines.txt
grep -oh "Otu[0-9]*" ${datadir}/number.lines.txt > ${datadir}/number.columns.txt

for nmbline in $(cat ${datadir}/number.columns.txt)
do
	grep -v "${nmbline}" ${datadir}/otu_table.txt > ${datadir}/otu_tmp.txt
	mv ${datadir}/otu_tmp.txt ${datadir}/otu_tmp2.txt
done

## Converter
biom convert \
	-i ${datadir}/otu_table.txt \
	-o ${datadir}/table.biom \
	--table-type="OTU table" \
	--to-json

## Importar no QIIME2
qiime tools import \
  --input-path ${datadir}/table.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV100Format \
  --output-path ${datadir}/feature-table.qza
  
## Gerar os valores de rarefação
echo "Summarizing counts..."
qiime feature-table summarize \
	    --i-table ${datadir}/feature-table.qza \
	    --o-visualization ${datadir}/feature-table.qzv \
	    --m-sample-metadata-file ${datadir}/metadata.txt
    
echo "Exporting HTML for counts..."
echo
## Exportando HTML
qiime tools export \
	--input-path ${datadir}/feature-table.qzv \
	--output-path ${datadir}/feature-table-html
	
## Pegar a rarefação
echo "Getting rarefaction value..."
rarevalue=`tail -1 ${datadir}/feature-table-html/sample-frequency-detail.csv | sed 's/.*,//' | sed 's/\..//'`
tail -15 ${datadir}/feature-table-html/sample-frequency-detail.csv
echo "Rarefaction is: ${rarevalue}"
echo

## Fazer a rarefação
qiime feature-table rarefy \
	--i-table ${datadir}/feature-table.qza \
	--p-sampling-depth ${rarevalue} \
	--o-rarefied-table ${datadir}/rare-table.qza
	
## Exportar o BIOM
echo "Exporting the rarefied BIOM..."
qiime tools export \
	--input-path ${datadir}/rare-table.qza \
	--output-path ${datadir}/rarefied
	
## Adicionando taxonomia
biom add-metadata \
	-i ${datadir}/rarefied/feature-table.biom \
	-o ${datadir}/rarefied/tabletax-rare.biom \
	--observation-metadata-fp ${datadir}/tax.tmp \
	--sc-separated taxonomy \
	--float-fields confidence
	
## Gerar por nível
summarize_taxa.py \
	--absolute_abundance \
	-i ${datadir}/rarefied/tabletax-rare.biom \
	-o ${datadir}/rarefied/tabulated
	
## Gerando a tabela de OTUs brutas
biom convert \
	-i ${datadir}/rarefied/tabletax-rare.biom \
	-o ${datadir}/rarefied/tabulated/tabletax-rare_L1_OTUs.parsed.txt \
	--to-tsv
	
for txtfile in $(ls ${datadir}/rarefied/tabulated/*.txt); do
	txtbase=$(basename ${txtfile} .txt)
	sed 's/.*;//g' ${txtfile} > ${datadir}/rarefied/tabulated/${txtbase}.parsed.txt
done

mv ${datadir}/rarefied/tabulated/*parsed* ${datadir}/output/
mv ${datadir}/* ./tmp/
mkdir -p ${processed}
mv ./tmp/* ${processed}/
mv ${processed}/output ${datadir}/
rm -rf ./tmp

rename "L2" "L2_phylum" ${datadir}/output/*
rename "L3" "L3_class" ${datadir}/output/*
rename "L4" "L4_order" ${datadir}/output/*
rename "L5" "L5_family" ${datadir}/output/*
rename "L6" "L6_genus" ${datadir}/output/*

## Copiar as amostras
cp ${processed}/feature-table-html/sample-frequency-detail.csv ${datadir}/output/samples.csv

## Re-cria as entradas
mkdir -p ${datadir}/input
cp ${processed}/metadata.txt ${datadir}/input/
sed 's/sample-id/#NAME/' ${datadir}/metadata.txt > ${datadir}/input/metadata.microbiomeAnalyst.txt
cp ${processed}/sintax.microbiomeAnalyst.txt ${datadir}/input/
cp ${processed}/otu_table.txt ${datadir}/input/

## Salva com o nome da análise
mkdir -p ./newtemp
mv ${datadir}/* ./newtemp
mkdir -p ${datadir}/${prepr_no}_${annot_no}_work_in_R
mv ./newtemp/* ${datadir}/${prepr_no}_${annot_no}_work_in_R
rm -rf ./newtemp

rename "tabletax-rare_" "" $datadir/${prepr_no}_${annot_no}_work_in_R/output/*
rename ".parsed" ""  $datadir/${prepr_no}_${annot_no}_work_in_R/output/*
 
echo "----------------------------------------------" >> ./analises/masterlog
echo "${merge_no}\tTentativa de juntar os pares" >> ./analises/masterlog
echo "${prepr_no}\tPré-processamento número" >> ./analises/masterlog
echo "${annot_no}\tTentativa de anotação número" >> ./analises/masterlog
echo "${annotation_parameter}\tOrganismo explorado" >> ./analises/masterlog
echo "${p_diffs}\tDiferença entre bases para juntar os pares" >> ./analises/masterlog
echo "${p_id}\tIdentidade para juntar os pares" >> ./analises/masterlog
echo "${p_ee}\tErros esperados" >> ./analises/masterlog
echo "${p_trunc}\tTamanho de truncagem" >> ./analises/masterlog
echo "${p_trimR1}\tTamanho de ponta da poda R1" >> ./analises/masterlog
echo "${p_trimR2}\tTamanho de ponta da poda R2" >> ./analises/masterlog
echo "----------------------------------------------" >> ./analises/masterlog


```

#### Adendo para ZOTUs

```shell
# datadir deve conter:
# otu_table.txt - matriz de OTUs
# sintax.microbiomeAnalyst.txt - arquivo de taxonomia formatado
# metadata.tsv - arquivo de metadados do QIIME2
source activate qiime2

## Copia a taxonomia
cp ./analises/E13_microbiome-analyst/${prepr_no}_ztaxonomy.txt ./analises/E14_work_in_R/sintax.microbiomeAnalyst.txt

## Copia a OTU-table
cp ./analises/E13_microbiome-analyst/${prepr_no}_zotu_table.txt ./analises/E14_work_in_R/otu_table.txt

## Copia os metadados
cp ./dados/intel/meta* ./analises/E14_work_in_R/metadata.txt

datadir="./analises/E14_work_in_R"
processed="${datadir}/processed"
output="${datadir}/output"
mkdir -p ${output}
mkdir -p ./tmp
 
## Teste: parse
sed 's/\t/;/g' ${datadir}/sintax.microbiomeAnalyst.txt | \
sed 's/\(Zotu[0-9]*\);/\1\t/g' | \
sed 's/$/\t1.00/g' |
sed 1d > ${datadir}/tax.tmp
sed -i '1 i\#OTUID\ttaxonomy\tconfidence' ${datadir}/tax.tmp

## Remoção de cloroplastos
grep -E "[cC]hloroplast" ${datadir}/tax.tmp > ${datadir}/chloro.lines.txt
grep -oh "Zotu[0-9]*" ${datadir}/chloro.lines.txt > ${datadir}/chloro.columns.txt

for chline in $(cat ${datadir}/chloro.columns.txt)
do
	grep -v "${chline}" ${datadir}/otu_table.txt > ${datadir}/otu_tmp.txt
	mv ${datadir}/otu_tmp.txt ${datadir}/otu_table.txt
done

## Remoção de mitocôndrias

grep -E "[mM]itochondria" ${datadir}/tax.tmp > ${datadir}/mito.lines.txt
grep -oh "Zotu[0-9]*" ${datadir}/mito.lines.txt > ${datadir}/mito.columns.txt

for mtline in $(cat ${datadir}/mito.columns.txt)
do
	grep -v "${mtline}" ${datadir}/otu_table.txt > ${datadir}/otu_tmp.txt
	mv ${datadir}/otu_tmp.txt ${datadir}/otu_table.txt
done

## Remoção de ruídos numéricos
grep -E "1\s1\s1\s1" ${datadir}/tax.tmp > ${datadir}/number.lines.txt
grep -E "2\s2\s2\s2" ${datadir}/tax.tmp > ${datadir}/number.lines.txt
grep -E "3\s3\s3\s3" ${datadir}/tax.tmp > ${datadir}/number.lines.txt
grep -oh "Zotu[0-9]*" ${datadir}/number.lines.txt > ${datadir}/number.columns.txt

for nmbline in $(cat ${datadir}/number.columns.txt)
do
	grep -v "${nmbline}" ${datadir}/otu_table.txt > ${datadir}/otu_tmp.txt
	mv ${datadir}/otu_tmp.txt ${datadir}/otu_tmp2.txt
done

## Converter
biom convert \
	-i ${datadir}/otu_table.txt \
	-o ${datadir}/table.biom \
	--table-type="OTU table" \
	--to-json

## Importar no QIIME2
qiime tools import \
  --input-path ${datadir}/table.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV100Format \
  --output-path ${datadir}/feature-table.qza
  
## Gerar os valores de rarefação
echo "Summarizing counts..."
qiime feature-table summarize \
	    --i-table ${datadir}/feature-table.qza \
	    --o-visualization ${datadir}/feature-table.qzv \
	    --m-sample-metadata-file ${datadir}/metadata.txt
    
echo "Exporting HTML for counts..."
echo
## Exportando HTML
qiime tools export \
	--input-path ${datadir}/feature-table.qzv \
	--output-path ${datadir}/feature-table-html
	
## Pegar a rarefação
echo "Getting rarefaction value..."
rarevalue=`tail -1 ${datadir}/feature-table-html/sample-frequency-detail.csv | sed 's/.*,//' | sed 's/\..//'`
tail -15 ${datadir}/feature-table-html/sample-frequency-detail.csv
echo "Rarefaction is: ${rarevalue}"
echo

## Fazer a rarefação
qiime feature-table rarefy \
	--i-table ${datadir}/feature-table.qza \
	--p-sampling-depth ${rarevalue} \
	--o-rarefied-table ${datadir}/rare-table.qza
	
## Exportar o BIOM
echo "Exporting the rarefied BIOM..."
qiime tools export \
	--input-path ${datadir}/rare-table.qza \
	--output-path ${datadir}/rarefied
	
## Adicionando taxonomia
biom add-metadata \
	-i ${datadir}/rarefied/feature-table.biom \
	-o ${datadir}/rarefied/tabletax-rare.biom \
	--observation-metadata-fp ${datadir}/tax.tmp \
	--sc-separated taxonomy \
	--float-fields confidence
	
## Gerar por nível
summarize_taxa.py \
	--absolute_abundance \
	-i ${datadir}/rarefied/tabletax-rare.biom \
	-o ${datadir}/rarefied/tabulated
	
## Gerando a tabela de OTUs brutas
biom convert \
	-i ${datadir}/rarefied/tabletax-rare.biom \
	-o ${datadir}/rarefied/tabulated/tabletax-rare_L1_OTUs.parsed.txt \
	--to-tsv
	
for txtfile in $(ls ${datadir}/rarefied/tabulated/*.txt); do
	txtbase=$(basename ${txtfile} .txt)
	sed 's/.*;//g' ${txtfile} > ${datadir}/rarefied/tabulated/${txtbase}.parsed.txt
done

mv ${datadir}/rarefied/tabulated/*parsed* ${datadir}/output/
mv ${datadir}/* ./tmp/
mkdir -p ${processed}
mv ./tmp/* ${processed}/
mv ${processed}/output ${datadir}/
rm -rf ./tmp

rename "L2" "L2_phylum" ${datadir}/output/*
rename "L3" "L3_class" ${datadir}/output/*
rename "L4" "L4_order" ${datadir}/output/*
rename "L5" "L5_family" ${datadir}/output/*
rename "L6" "L6_genus" ${datadir}/output/*

## Copiar as amostras
cp ${processed}/feature-table-html/sample-frequency-detail.csv ${datadir}/output/samples.csv

## Re-cria as entradas
mkdir -p ${datadir}/input
cp ${processed}/metadata.txt ${datadir}/input/
sed 's/sample-id/#NAME/' ${datadir}/metadata.txt > ${datadir}/input/metadata.microbiomeAnalyst.txt
cp ${processed}/sintax.microbiomeAnalyst.txt ${datadir}/input/
cp ${processed}/otu_table.txt ${datadir}/input/

## Salva com o nome da análise
mkdir -p ./newtemp
mv ${datadir}/* ./newtemp
mkdir -p ${datadir}/${prepr_no}_${annot_no}_work_in_R
mv ./newtemp/* ${datadir}/${prepr_no}_${annot_no}_work_in_R
rm -rf ./newtemp

rename "tabletax-rare_" "" $datadir/${prepr_no}_${annot_no}_work_in_R/output/*
rename ".parsed" ""  $datadir/${prepr_no}_${annot_no}_work_in_R/output/*
 
echo "----------------------------------------------" >> ./analises/masterlog
echo "${merge_no}\tTentativa de juntar os pares" >> ./analises/masterlog
echo "${prepr_no}\tPré-processamento número" >> ./analises/masterlog
echo "${annot_no}\tTentativa de anotação número" >> ./analises/masterlog
echo "${annotation_parameter}\tOrganismo explorado" >> ./analises/masterlog
echo "${p_diffs}\tDiferença entre bases para juntar os pares" >> ./analises/masterlog
echo "${p_id}\tIdentidade para juntar os pares" >> ./analises/masterlog
echo "${p_ee}\tErros esperados" >> ./analises/masterlog
echo "${p_trunc}\tTamanho de truncagem" >> ./analises/masterlog
echo "${p_trimR1}\tTamanho de ponta da poda R1" >> ./analises/masterlog
echo "${p_trimR2}\tTamanho de ponta da poda R2" >> ./analises/masterlog
echo "----------------------------------------------" >> ./analises/masterlog
```

## Etapa adicional (15): Anotação com PICRUSt

Experimental

```shell
## Implementado como default-picrust

picrust2_pipeline.py \
	-s ./analises/E08_pick_otus_zotus/otus.P01.fa \
	-i ./analises/E09_otu_tables/otutab.P01.tsv \
	-o ./analises/E15_PICRUSt2 \
	-p 20
	
## Com parâmetros
picrust2_pipeline.py \
	-s ./analises/E08_pick_otus_zotus/otus.${prepr_no}.fa \
	-i ./analises/E09_otu_tables/otutab.${prepr_no}.tsv \
	-o ./analises/E15_PICRUSt2.${prepr_no} \
	-p 20
```

## E16. Diversidade-alfa

```shell
usearch11 \
	-alpha_div ./analises/E09_otu_tables/otutab.${prepr_no}.tsv \
	-output ./analises/E16_alpha/alpha01.default.txt
	
usearch11 \
	-alpha_div ./analises/E09_otu_tables/otutab.${prepr_no}.tsv \
	-output ./analises/E16_alpha/alpha02.gini.txt \
	-metrics gini_simpson
	
usearch11 \
	-alpha_div ./analises/E09_otu_tables/otutab.${prepr_no}.tsv \
	-output ./analises/E16_alpha/alpha03.chao1_berger.txt
```



## E17. Diversidade-beta



## Adendo III: Salvando os parâmetros da corrida

```shell
echo "----------------------------------------------" >> ./analises/masterlog
echo "${merge_no}\tTentativa de juntar os pares" >> ./analises/masterlog
echo "${prepr_no}\tPré-processamento número" >> ./analises/masterlog
echo "${annot_no}\tTentativa de anotação número" >> ./analises/masterlog
echo "${p_organism}\tOrganismo explorado" >> ./analises/masterlog
echo "${p_diffs}\tDiferença entre bases para juntar os pares" >> ./analises/masterlog
echo "${p_id}\tIdentidade para juntar os pares" >> ./analises/masterlog
echo "${p_ee}\tErros esperados" >> ./analises/masterlog
echo "${p_trunc}\tTamanho de truncagem" >> ./analises/masterlog
echo "${p_trimR1}\tTamanho de ponta da poda R1" >> ./analises/masterlog
echo "${p_trimR2}\tTamanho de ponta da poda R2" >> ./analises/masterlog
echo "----------------------------------------------" >> ./analises/masterlog

```

## Adendo IV: Tamanho das regiões 16S

| Região sequenciada | Tamanho (pb) |
| ------------------ | ------------ |
| V2-V3              | 205          |
| V3-V4              | 415          |

## Relatórios de qualidade

```shell
## Diretório de logs
mkdir -p reports

## Salva o número de sequencias
grep -c "^@M" dados/brutos/*R1* >> reports/fq.log.1c.txt

for sample in $(cut -f1 ./dados/intel/metadados.tsv | tac | head -n -1 | tac);
do
	echo "${sample}"
	echo
	echo
	
	grep -c "^@${sample}" analises/E02_mergepairs/merged.M01.fq >> reports/merged.log.1c.txt
	grep -c "^@${sample}" analises/E03_trimming/trimmed.M01.${p_trunc}.fq >> reports/trimmed.log.1c.txt
	grep -c "^>${sample}" analises/E05_filter/filtered.P01.EE_${p_ee}.fa >> reports/filtered.log.1c.txt
	grep -c "^${sample}" analises/E09_otu_tables/mapout.P01.tsv	>> reports/mapped.log.1c.txt
done


echo "GG"
echo


for sample in $(cut -f1 ./dados/intel/metadados.tsv | tac | head -n -1 | tac)
do
	# grep -c "^@${sample}" analises/E02_mergepairs/merged.M01.fq
	# grep -c "^@${sample}" analises/E03_trimming/trimmed.M01.${p_trunc}.fq
	# grep -c "^>${sample}" analises/E05_filter/filtered.P01.EE_${p_ee}.fa
	grep -c "^${sample}" analises/E09_otu_tables/mapout.P01.tsv
done
```
