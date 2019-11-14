## Rodando DESeq2 time series analysis
# Autora: Maria Fernanda Zaneli Campanari (nandacampanari@gmail.com)

# Declarando o caminho
setwd("/usr/local/data/genome/fernandacampanari/giovanna/dados")

Matrix.genes <- read.delim(file = "/usr/local/data/genome/fernandacampanari/giovanna/dados/counts_int.txt")
length(Matrix.genes)
Matrix.genes

# Lendo a matrix e cnvertando os números para inteiros
Matrix.genes <- read.delim("/usr/local/data/genome/fernandacampanari/giovanna/dados/counts_int.txt", row.names=1)
Matrix.genes <- round(Matrix.genes) 

# Modificando ordem das colunas
Matrix.genes[,c(13,14,15,4,5,6,7,8,9,10,11,12,1,2,3)]
View(Matrix.genes[,c(13,14,15,4,5,6,7,8,9,10,11,12,1,2,3)])

# Guardando a matrix de genes corrigida (magecor)
magecor <- Matrix.genes[,c(13,14,15,4,5,6,7,8,9,10,11,12,1,2,3)]

# Chamando tabela com ordem temporal
Id_times <- read.delim(
                      "/usr/local/data/genome/fernandacampanari/giovanna/dados/samples_deseq.txt",
                      stringsAsFactors = F,
                      row.names=1)

#ctsgene <- as.matrix(Matrix.genes, sep ="\t", row.names="gene_id")

# Modificando o nome da tabela de ordem temporal
coldata <- Id_times

# Modificando a coluna de id
#row.names(ctsgene) <- ctsgene[,1]
#ctsgene <-ctsgene[,-1]

## Rodando o Deseq
# bilioteca DESeq

library("DESeq2")

# atribuindo o nome de dds (DESeqdataset)
dds <- DESeqDataSetFromMatrix(countData = magecor,
                              colData = coldata,
                              design = ~ time)
dds

# Fazendo o teste Likelihood ratio test

dds <- DESeq(dds, test="LRT", reduced=~1)

#chamando os resultados
res <- results(dds)

#chamando os resultados com pajustado
res.padj <- subset(res, padj < 0.05)

# Guardando no servidor

write.table(as.matrix(res.padj),
            file="/usr/local/data/genome/fernandacampanari/giovanna/analises/2019-10-25_deseq/deseq.tsv",
            quote = F,
            sep = '\t',
            row.names = F)

View(as.matrix(res.padj))

res.padj.t5vt1 <- subset(res.padj, log2FoldChange > 2)

# Exportar a matriz de contagens normalizadas
dds <- estimateSizeFactors(dds)
deseq.counts = counts(dds, normalized=TRUE)

# Extrair os DEs
deseq.counts.padj = subset(deseq.counts,
                           row.names(deseq.counts) %in% row.names(res.padj))

# Exportando para o servidor
write.table(as.matrix(deseq.counts.padj),
            file="/usr/local/data/genome/fernandacampanari/giovanna/analises/2019-10-25_deseq/deseq.counts.tsv",
            quote = F,
            sep = '\t',
            row.names = T)
