# Script modelo para testes de hipóteses paramétrico e não-paramétrico

#### [.. 1 Ativação dos pacotes ..] ####
install.packages(c(
  "emmeans", "readxl", "caret",
  "lmtest", "car", "rcompanion",
  "userfriendlyscience", "pgirmess"
))

library(emmeans)
library(readxl)
library(caret)
library(lmtest)
library(car)
library(rcompanion)
library(userfriendlyscience)
library(pgirmess)
library(multcomp)
library(multcompView)
library(ggplot2)
library(tidyverse)
library(gapminder)
library(forcats)

#### [.. 2. Validação dos dados ..] ####

tabela <- read.table(file = "clipboard", sep = "\t", header = T)
tabela.bkp <- tabela

# Distribuições
# 1. Parece binomial
head(tabela)

# 2. Parece normal, mas não é
hist(tabela$r)
hist(tabela$pa)

#### 3. Testes de Shapiro ####
# Não vem de normal!
shapiro.raiz <- shapiro.test(sort(tabela$r))
shapiro.aere <- shapiro.test(sort(tabela$pa))

# Rejeição no teste de Shapiro:
# Os dados não vem de modelo normal

#### 4. Anova forçada ####
anova.fact01.r <- Anova(aov(r~trat + bloco, tabela))
anova.fact02.pa <- Anova(aov(pa~trat + bloco, tabela))

fact01.r.dunnett <- cld(emmeans(aov(r~trat, tabela),
                                specs=c("trat"), adjust="dunnettx"),
                        alpha=0.05, Letters=letters, adjust="dunnettx")

fact02.pa.dunnett <- cld(emmeans(aov(pa~trat, tabela),
                                specs=c("trat"), adjust="dunnettx"),
                        alpha=0.05, Letters=letters, adjust="dunnettx")


#### 5. Kruskal-Wallis ####
kruskal.test(formula = r~trat, data = tabela)
kruskal.test(formula = r~bloco, data = tabela)
kruskal.test(formula = r~bact, data = tabela)
kruskal.test(formula = r~metal, data = tabela)
kruskal.test(formula = r~bact.metal, data = tabela)
kruskal.test(formula = r~which.metal, data = tabela)


kruskal.test(formula = pa~trat, data = tabela)
kruskal.test(formula = pa~bloco, data = tabela)
kruskal.test(formula = pa~bact, data = tabela)
kruskal.test(formula = pa~metal, data = tabela)
kruskal.test(formula = pa~bact.metal, data = tabela)
kruskal.test(formula = pa~which.metal, data = tabela)

# estrinchamento...
kruskal.r <- kruskalmc(r~trat, data = tabela)
kruskal.pa <- kruskalmc(pa~trat, data = tabela)
kruskal.me <- kruskalmc(pa~metal, data = tabela)

# Bactéria
kruskal.r.bact <- kruskalmc(r~bact, data = tabela)
kruskal.pa.bact <- kruskalmc(pa~bact, data = tabela)

# Metal
kruskal.r.metal <- kruskalmc(r~metal, data = tabela)
kruskal.pa.metal <- kruskalmc(pa~metal, data = tabela)

# Bactéria + metal
kruskal.r.bact.metal <- kruskalmc(r~bact.metal, data = tabela)
kruskal.pa.bact.metal <- kruskalmc(pa~bact.metal, data = tabela)

# Qual metal
kruskal.r.which.metal <- kruskalmc(r ~ which.metal, data = tabela)
kruskal.pa.which.metal <- kruskalmc(pa ~ which.metal, data = tabela)

#### 6. Games-Howell (teste) ####

ph.fact01.r <- oneway(y = tabela$r, x = tabela$trat,
                      posthoc = 'games-howell', levene = T,
                      posthocLetters = T, pvalueDigits = 4,
                      fullDescribe = T, plot = T)

#### 7. Figuras ####

## Por metal, raízes
ggplot(tabela, aes(x = fct_reorder(trat, r), y = r, fill = trat)) +  geom_boxplot() + geom_jitter(width=0.1, alpha=0.2)

ggplot(tabela, aes(x = fct_reorder(trat, pa), y = pa, fill = trat)) +  geom_boxplot() + geom_jitter(width=0.1, alpha=0.2)

# Somente bactéria
## Raízes
ggplot(tabela, aes(x = fct_reorder(bact, r), y = r, fill = bact)) +  geom_boxplot() + geom_jitter(width=0.1, alpha=0.2)

## Parte aérea
ggplot(tabela, aes(x = fct_reorder(bact, pa), y = pa, fill = bact)) + geom_boxplot() + geom_jitter(width=0.1, alpha=0.2)

# Somente metal
ggplot(tabela, aes(x = fct_reorder(metal, r),
                   y = r,
                   fill = metal)) + 
  geom_boxplot() + geom_jitter(width=0.1, alpha=0.2)

# Bactéria e metal r
ggplot(tabela, aes(x = fct_reorder(bact.metal, r),
                   y = r,
                   fill = bact.metal)) + 
  geom_boxplot() + geom_jitter(width=0.1, alpha=0.2)

# Bactéria e metal pa
ggplot(tabela, aes(x = fct_reorder(bact.metal, pa),
                   y = pa,
                   fill = bact.metal)) + 
  geom_boxplot() + geom_jitter(width=0.1, alpha=0.2)

# Qual metal?
ggplot(tabela, aes(x = fct_reorder(which.metal, pa),
                   y = pa, fill = which.metal)) + 
  geom_boxplot() + geom_jitter(width=0.1, alpha=0.2) +
  xlab("Metais") +
  ylab("Comprimento da parte aérea") +
  ggtitle("Parte aérea em resposta a qual metal")

ggplot(tabela, aes(x = fct_reorder(which.metal, r),
                   y = r, fill = which.metal)) + 
  geom_boxplot() + geom_jitter(width=0.1, alpha=0.2) +
  xlab("Metais") +
  ylab("Comprimento da raiz") +
  ggtitle("Parte radicular em resposta a qual metal")

#### Gerar as letras ####

agg.sum <- aggregate(tabela[, c(2:3)], list(tabela$trat), mean)
agg.metal <- aggregate(tabela[, c(1:2)], list(tabela$which.metal), mean)

# agg.sum <- agg.sum[,-2]
colnames(agg.sum)[1] = "trat"

agg.r <- agg.sum[order(agg.sum$r),]
agg.pa <- agg.sum[order(agg.sum$pa),]

## Reordenar fator
reorder(tabela$trat, by = as.factor(agg.r$Group.1))

#### 8. Post-hoc para Kruskal-Wallis: Mann Withey U ####
PT.r <- pairwise.wilcox.test(tabela$r, tabela$trat)
PT.pa <- pairwise.wilcox.test(tabela$pa, tabela$trat)

## Post-hoc para qual metal ##
PT.r.metal <- pairwise.wilcox.test(tabela$r, tabela$which.metal)
PT.pa.metal <- pairwise.wilcox.test(tabela$pa, tabela$which.metal)

PT.r = PT.r$p.value
PT.pa = PT.pa$p.value
PT1.r = fullPTable(PT.r)
PT1.pa = fullPTable(PT.pa)

### Letras para raízes
r.letters <- multcompLetters(PT1.r, compare = "<", threshold = 0.05, Letters = letters, reversed = F)
pa.letters <- multcompLetters(PT1.pa, compare = "<", threshold = 0.05, Letters = letters, reversed = F)

r.WM.letters <- multcompLetters(PT.r.metal, compare = "<", threshold = 0.05, Letters = letters, reversed = F)

## Cria a matriz de letras
df.r <- data.frame("trat" = names(r.letters$Letters), "r" = r.letters$Letters)
colnames(df.r)[2] = "posthoc.r"
tabela <- merge(tabela, df.r, by = "trat")

df.pa <- data.frame("trat" = names(pa.letters$Letters), "pa" = pa.letters$Letters)
colnames(df.pa)[2] = "posthoc.pa"
tabela <- merge(tabela, df.pa, by = "trat")

## Tabela de significativos
agg.table <- merge(agg.sum, df.r, by = "trat")
agg.table <- merge(agg.table, df.pa, by = "trat")

#### 9. Testes no Cromo ####
kruskal.test(formula = r~crome.test, data = tabela)
kruskal.test(formula = pa~crome.test, data = tabela)

kruskal.crome.ra <- kruskalmc(r~crome.test, data = tabela)
kruskal.crome.pa <- kruskalmc(pa~crome.test, data = tabela)

ggplot(tabela, aes(x = fct_reorder(crome.test, pa),
                   y = pa, fill = crome.test)) + 
  geom_boxplot() + geom_jitter(width=0.1, alpha=0.2) +
  xlab("Metais") +
  ylab("Comprimento da parte aérea") +
  ggtitle("Cromo em relação a concentrações, e presença ou não de bactéria")
