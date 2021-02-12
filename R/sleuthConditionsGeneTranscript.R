# 2020-02-08 First edit
# Sleuth: Transcriptoma
# Rodar do diretório "04_kallisto"

# Bibliotecas
library("sleuth")
library("dplyr")

## Escolhe as duas condições a serem rodadas
condition_one = "MI-J4"
condition_two = "MI-J3"

# Diretório
base_dir <- "."

# Nomes das amostras a serem analisadas

sample_id = c("SRR5684403","SRR5684404","SRR5684405","SRR5684406","SRR5684407","SRR5684408","SRR5684409","SRR5684410","SRR5684411","SRR5684412","SRR5684413","SRR5684414","SRR5684415","SRR5684416","SRR5684417")

# Lista dos diretórios (paste)
# lista dos diretorios criados pelo kallisto (output). Lembrar de manter a msm ordem q no sample_id

## Gerar todos de uma vez
kallisto_out = c(
	file.path(base_dir, "SRR5684403_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684404_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684405_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684406_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684407_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684408_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684409_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684410_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684411_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684412_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684413_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684414_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684415_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684416_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684417_kallisto_transcripts_index.out")
)

# Cola os nomes da amostra nos diretórios do Kallisto
names(kallisto_out) <- sample_id

# Ler o desenho experimental
s2c <- read.table(file.path(base_dir, "matrix_kallisto", "design_matrix.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::mutate(s2c, path = kallisto_out)

# Adiciona o arquivo "transcrito para gene", t2gene
## Com "fill" true: preenche os transcritos vazios
t2g <- read.table(file.path(base_dir, "t2gene.txt"), header = TRUE, fill = TRUE, stringsAsFactors=FALSE)

# Sleuth propriamente dito
## Calculo (da condição), normaliza
so.transcripts <- sleuth_prep(s2c, ~condition, target_mapping = t2g, extra_bootstrap_summary = TRUE)

## Calcula sobre os genes
so.gene <- sleuth_prep(s2c, ~condition, target_mapping = t2g, aggregation_column= "gene",
     extra_bootstrap_summary = TRUE, max_bootstrap = 10, normalize = TRUE, gene_mode = TRUE)

## Cálculo da expressão diferencial
so.transcripts <- sleuth_fit(so.transcripts, num_cores=10)
so.gene <- sleuth_fit(so.gene)

## Teste de wald
so.transcripts <- sleuth_wt(so.transcripts, "condition")
models(so.transcripts)

## Grava os resultados
results.so.tr <- sleuth_results(so.transcripts, test="conditionMI-Female", test_type="wald")

## Filtra os diferenciais (transcrito)
sl.sig.tr <- dplyr::filter(
    results.so.tr,
    qval <= 0.05
)

# Output table
## Apenas um exemplo, troque para as outras condições
write.table(
	x = as.data.frame(sl.sig.tr),
	file = "./out.MI.female.tsv",
	quote = FALSE,
	sep = "\t",
	row.names = FALSE
)



























# 2020-02-12 Testes
# Bibliotecas
library("sleuth")
library("dplyr")

## Escolhe as duas condições a serem rodadas
condition_one = "MI-J4"
condition_two = "MI-J3"

# Diretórios
setwd("/usr/local/data/genome/fernandacampanari/area_gi/ZZ_main/analysis/04_kallisto")
base_dir <- "."

# Nomes das amostras a serem analisadas

sample_id = c("SRR5684403","SRR5684404","SRR5684405","SRR5684406","SRR5684407","SRR5684408","SRR5684409","SRR5684410","SRR5684411","SRR5684412","SRR5684413","SRR5684414","SRR5684415","SRR5684416","SRR5684417")

# Lista dos diretórios (paste)
# lista dos diretorios criados pelo kallisto (output). Lembrar de manter a msm ordem q no sample_id

## Gerar todos de uma vez
kallisto_out = c(
	file.path(base_dir, "SRR5684403_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684404_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684405_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684406_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684407_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684408_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684409_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684410_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684411_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684412_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684413_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684414_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684415_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684416_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684417_kallisto_transcripts_index.out")
)

# Cola os nomes da amostra nos diretórios do Kallisto
names(kallisto_out) <- sample_id

# Ler o desenho experimental
s2c <- read.table(file.path(base_dir, "matrix_kallisto", "design_matrix.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::mutate(s2c, path = kallisto_out)

# Adiciona o arquivo "transcrito para gene", t2gene
## Com "fill" true: preenche os transcritos vazios
t2g <- read.table(file.path(base_dir, "t2gene.txt"), header = TRUE, fill = TRUE, stringsAsFactors=FALSE)

## Renivela as condições
s2c$condition <- as.factor(s2c$condition)

#### [. PARTE NOVA . ] ####
## Faz o nivelamento e renivelamento
### Escolhe uma condição de base, deve ser fixado, tudo será comparado com ele
### Roda somente para essa condição
s2c$condition <- as.factor(s2c$condition)
s2c$condition <- relevel(s2c$condition, condition_one)

### Sleuth propriamente dito
so.transcripts <- sleuth_prep(s2c, ~condition, target_mapping = t2g, extra_bootstrap_summary = TRUE, num_cores=10)
so.transcripts <- sleuth_fit(so.transcripts, num_cores=10)

# Caso queira checar as opções de resultados
models(so.transcripts)

### Escolhe a condição pareada
paired_condition = paste0("condition", condition_two)
so.transcripts <- sleuth_wt(so.transcripts, paired_condition)
so.results <- sleuth_results(so.transcripts, test=paired_condition, test_type="wald")

### Gera os significativos
diff.sig <- dplyr::filter(so.results, qval <= 0.05)

### Gera o nome do arquivo de saida
output_name <- paste0("./out.", condition_one, "_vs_", condition_two, ".tsv")

write.table(
	x = as.data.frame(diff.sig),
	file = output_name,
	quote = FALSE,
	sep = "\t",
	row.names = FALSE
)






































#### Função para genes
# 2020-02-12 Testes
# Bibliotecas
library("sleuth")
library("dplyr")

# Diretórios
setwd("/usr/local/data/genome/fernandacampanari/area_gi/ZZ_main/analysis/04_kallisto")
base_dir <- "."

## Escolhe as duas condições a serem rodadas
condition_one = "MI-J4"
condition_two = "MI-Female"

# Nomes das amostras a serem analisadas

sample_id = c("SRR5684403","SRR5684404","SRR5684405","SRR5684406","SRR5684407","SRR5684408","SRR5684409","SRR5684410","SRR5684411","SRR5684412","SRR5684413","SRR5684414","SRR5684415","SRR5684416","SRR5684417")

# Lista dos diretórios (paste)
# lista dos diretorios criados pelo kallisto (output). Lembrar de manter a msm ordem q no sample_id

## Gerar todos de uma vez
kallisto_out = c(
	file.path(base_dir, "SRR5684403_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684404_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684405_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684406_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684407_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684408_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684409_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684410_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684411_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684412_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684413_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684414_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684415_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684416_kallisto_transcripts_index.out"),
	file.path(base_dir, "SRR5684417_kallisto_transcripts_index.out")
)

# Cola os nomes da amostra nos diretórios do Kallisto
names(kallisto_out) <- sample_id

# Ler o desenho experimental
s2c <- read.table(file.path(base_dir, "matrix_kallisto", "design_matrix.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::mutate(s2c, path = kallisto_out)

# Adiciona o arquivo "transcrito para gene", t2gene
## Com "fill" true: preenche os transcritos vazios
t2g <- read.table(file.path(base_dir, "t2gene.txt"), header = TRUE, fill = TRUE, stringsAsFactors=FALSE)

## Renivela as condições
s2c$condition <- as.factor(s2c$condition)

### Roda somente para essa condição
s2c$condition <- as.factor(s2c$condition)
s2c$condition <- relevel(s2c$condition, condition_one)

### Sleuth propriamente dito
## Calcula sobre os genes
so.gene <- sleuth_prep(s2c, ~condition, target_mapping = t2g, aggregation_column= "gene",
     extra_bootstrap_summary = TRUE, max_bootstrap = 10, normalize = TRUE, gene_mode = TRUE)

so.gene <- sleuth_fit(so.gene, num_cores=10)

# Caso queira checar as opções de resultados
models(so.gene)

### Escolhe a condição pareada
#### Escolher o tratamento que irá VS o base
#### Depois, deve ser trocado pra J3, Egg, ...
paired_condition = paste0("condition", condition_two)
so.gene <- sleuth_wt(so.gene, paired_condition)
so.results <- sleuth_results(so.gene, test=paired_condition, test_type="wald")

### Gera os significativos
diff.sig <- dplyr::filter(so.results, qval <= 0.05)

### Gera o nome do arquivo de saida
output_name <- paste0("./gene.out.", condition_one, "_vs_", condition_two, ".tsv")

write.table(
	x = as.data.frame(diff.sig),
	file = output_name,
	quote = FALSE,
	sep = "\t",
	row.names = FALSE
)






