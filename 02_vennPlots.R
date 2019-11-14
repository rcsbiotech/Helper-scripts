# Diagramas de Venn em diferentes estilos/número de conjuntos

install.packages('VennDiagram')
install.packages('venn')
library(venn)
library(VennDiagram)

## Teste: venn
set.seed(12345)
x <- list(
  First = 1:20, 
  Second = 10:30, 
  Third = sample(25:50, 15), 
  Fourth = 12:15, 
  Fifth = 3)

venn(x)

## Passa a pasta
setwd("/work/rcsilva/projects/helper/2019-02-25_cnptia/data/de")

#### Etapa 1. Importando todos os dados como lista ####
gls <- read.table("gls.txt", stringsAsFactors = F)
int <- read.table("int.txt", stringsAsFactors = F)
lrv <- read.table("lrv.txt", stringsAsFactors = F)
nnf <- read.table("nnf.txt", stringsAsFactors = F)
ovr <- read.table("ovr.txt", stringsAsFactors = F)

### Concatena todas as isoformas em uma única matriz
all <- unique(rbind(gls, int, lrv, nnf, ovr))

venn.df <- data.frame(
  "Isoforms" = all,
  "GLS" = NA,
  "INT" = NA,
  "LRV" = NA,
  "NNF" = NA,
  "OVR" = NA
)

names(venn.df)[1] <- "Isoform"
names(venn.df) <- c("Isoform", "GLS (137)", "INT (134)", "LRV (14346)", "NNF (1051)", "OVR (251)")

#### Matriz VF ####
counter <- 1

for (isoform in venn.df$Isoform) {
  
  # GLS
  if (isoform %in% gls$V1) {
    venn.df$GLS[counter] <- 1
  } else {
    venn.df$GLS[counter] <- 0
  }
  
  #INT
  if (isoform %in% int$V1) {
    venn.df$INT[counter] <- 1
  } else {
    venn.df$INT[counter] <- 0
  }
  
  # LRV
  if (isoform %in% lrv$V1) {
    venn.df$LRV[counter] <- 1
  } else {
    venn.df$LRV[counter] <- 0
  }
  
  # NNF 
  if (isoform %in% nnf$V1) {
    venn.df$NNF[counter] <- 1
  } else {
    venn.df$NNF[counter] <- 0
  }
  
  # OVR
  if (isoform %in% ovr$V1) {
    venn.df$OVR[counter] <- 1
  } else {
    venn.df$OVR[counter] <- 0
  }
  
  counter <- counter + 1
}


## Cria o venn dos overlaps para DE
Figure <- venn(venn.df[,2:6],
               size=5000,
               ellipse=T,
               borders=F,
               zcolor="style")

plotname <- "/work/rcsilva/projects/helper/2019-02-25_cnptia/fig1.jpg"
ggsave(filename=plotname,
       plot=Figure)


#### 2. DE somente para GLS, OVR e INT ####
afew <- unique(rbind(gls, int, ovr))

trio.df <- data.frame(
  "Isoforms" = afew,
  "GLS" = NA,
  "INT" = NA,
  "OVR" = NA
)

names(trio.df)[1] <- "Isoform"

#### Matriz VF ####
counter <- 1

for (isoform in trio.df$Isoform) {
  
  # GLS
  if (isoform %in% gls$V1) {
    trio.df$GLS[counter] <- 1
  } else {
    trio.df$GLS[counter] <- 0
  }
  
  #INT
  if (isoform %in% int$V1) {
    trio.df$INT[counter] <- 1
  } else {
    trio.df$INT[counter] <- 0
  }
  
  # OVR
  if (isoform %in% ovr$V1) {
    trio.df$OVR[counter] <- 1
  } else {
    trio.df$OVR[counter] <- 0
  }
  
  counter <- counter + 1
}

## Figure
venn(
  trio.df[,2:4],
  borders=F,
  zcolor="style",
  ellipse=F
)
  
#### Etapa 3. Para os transcritos (notDE) ####

## Passa a pasta
setwd("/work/rcsilva/projects/helper/2019-02-25_cnptia/data/trans")

#### Etapa 1. Importando todos os dados como lista ####
gls <- read.table("gls.txt", stringsAsFactors = F)
int <- read.table("int.txt", stringsAsFactors = F)
lrv <- read.table("lrv.txt", stringsAsFactors = F)
nnf <- read.table("nnf.txt", stringsAsFactors = F)
ovr <- read.table("ovr.txt", stringsAsFactors = F)

### Concatena todas as isoformas em uma única matriz
all <- unique(rbind(gls, int, lrv, nnf, ovr))

venn.df <- data.frame(
  "Isoforms" = all,
  "GLS" = NA,
  "INT" = NA,
  "LRV" = NA,
  "NNF" = NA,
  "OVR" = NA
)

names(venn.df)[1] <- "Isoform"
names(venn.df) <- c("Isoform", "GLS", "INT", "LRV", "NNF", "OVR")

#### Código nesse formato muito lerdo - pulado ####
counter <- 1

for (isoform in venn.df$Isoform) {
  
  # GLS
  if (isoform %in% gls$V1) {
    venn.df$GLS[counter] <- 1
  } else {
    venn.df$GLS[counter] <- 0
  }
  
  #INT
  if (isoform %in% int$V1) {
    venn.df$INT[counter] <- 1
  } else {
    venn.df$INT[counter] <- 0
  }
  
  # LRV
  if (isoform %in% lrv$V1) {
    venn.df$LRV[counter] <- 1
  } else {
    venn.df$LRV[counter] <- 0
  }
  
  # NNF 
  if (isoform %in% nnf$V1) {
    venn.df$NNF[counter] <- 1
  } else {
    venn.df$NNF[counter] <- 0
  }
  
  # OVR
  if (isoform %in% ovr$V1) {
    venn.df$OVR[counter] <- 1
  } else {
    venn.df$OVR[counter] <- 0
  }
  
  counter <- counter + 1
}

#### Vetorizado ####
venn.df$GLS <- venn.df$Isoform %in% gls$V1
venn.df$INT <- venn.df$Isoform %in% int$V1
venn.df$LRV <- venn.df$Isoform %in% lrv$V1
venn.df$NNF <- venn.df$Isoform %in% nnf$V1
venn.df$OVR <- venn.df$Isoform %in% ovr$V1

## Cria o venn dos overlaps para DE
Figure <- venn(venn.df[,2:6],
               size=5000,
               ellipse=T,
               borders=F,
               zcolor="style")

#### DE somente para GLS, OVR e INT ####
afew <- unique(rbind(gls, int, ovr))

trio.df <- data.frame(
  "Isoforms" = afew,
  "GLS" = NA,
  "INT" = NA,
  "OVR" = NA
)

names(trio.df)[1] <- "Isoform"

trio.df$GLS <- trio.df$Isoform %in% gls$V1
trio.df$INT <- trio.df$Isoform %in% int$V1
trio.df$OVR <- trio.df$Isoform %in% ovr$V1

## Figure
venn(
  trio.df[,2:4],
  borders=F,
  zcolor="style",
  ellipse=F
)


