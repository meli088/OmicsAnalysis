getwd()
install.packages("gplots")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
BiocManager::install("edgeR")
a

library(limma)
library(edgeR)
library(gplots)

# Exploration des données
SamplesInfo <- read.table("TP2_Omiques_InfoSamples.tsv", h=T, sep="\t", stringsAsFactors = F)
CountData <- read.delim("TP2_Omiques_TableauComptage.tsv", row.names = 1)

head(CountData)
summary(CountData)

# Affichez les premières lignes du tableau de comptage, puis faites un graphique pour afficher la distribution du comptage de reads pour un des échantillons.
hist(log2(CountData[,1] +1 ), main = 'Histogramme dun ech', breaks = 50, xlim=c(-2, max(log2(CountData[,1] + 1))))
boxplot(log2(CountData+1), las=2, main='boxplot de tt')

# Trait c'est la médiane : à 0
# Une cellule n'exprime pas ts les gènes, ici y'a ttes les features de ts les gènes exprimés dans la cellule
# Bcp de gènes à 0
# 50% des données dans le IQR, partir de Q3, + 1,5*IQR : valeur seuil qui correspond à l'attendu th des 95%
# Grosses bandes noires : ya bcp de gènes exprimés par là mais on s'en blc

str(SamplesInfo)
# Convertir les colonnes en facteurs
SamplesInfo$condition <- factor(SamplesInfo$condition, levels=c("basal.virgin", "basal.pregnant", "basal.lactate", "luminal.virgin", "luminal.pregnant", "luminal.lactate")) # 6 conditions au total
SamplesInfo$cellType <- factor(SamplesInfo$cellType, levels=c("basal", "luminal"))
SamplesInfo$mother <- factor(SamplesInfo$mother, levels=c("virgin", "pregnant", "lactate"))
str(SamplesInfo)

# Convertir les 0 en NA, les supprimer
for(i in 1:ncol(CountData)){
  CountData[,i][which(CountData[,i]==0)] <- NA # retourne index == 0
}
CountData <- na.omit(CountData)
nrow(CountData)

# Créer un dataset utilisable pour la suite du TP (dgeFull) et un tableau de comptage transformé pour la visualisation.
# Création d'un tableau de comptage edgeR (objet DGElist)
dgeFull <- edgeR:: DGEList(CountData, remove.zeros = TRUE)

# Annotation du DGElist
dgeFull$samples$condition <- SamplesInfo$condition
dgeFull$samples$cellType <- SamplesInfo$cellType
dgeFull$samples$Mother <- SamplesInfo$mother

# Tableau de pseudo-comptage
pseudoCounts <- log10(dgeFull$counts + 1)

# Histogramme : entre 100 et 1000 lectures par gènes
hist(pseudoCounts[ ,"MCL1.DG"], main = "MCL1.DG", xlab = "nombre de reads par gène (échelle log10)", ylab="nombre de gènes")

# Boxplot : 2 réplicats mais comme même diff, préserver variabilité intra ech mais on veut que globalement médiane au même endroit et variance égale (boîte même taille)
boxplot(pseudoCounts, col = "red", las = 3, cex.names = 1, ylab="nombre de reads par gène (échelle log10)")
abline(h=median(pseudoCounts), col='blue', lwd=2, lty='dashed')
# La normalisation on y est pas encore car ttes les médianes ne sont pas alignées...
# Malgrè le fait que ce soit pas normalisé on s'att à ce qu'entre 2 duplicats ce soit la même...z
# Donc on va faire une MDS : calculer expr moyenne de chaq gène de chaq ech et regroupé les ind en fonction ressemblance
# MDS : ds cet exemple, distance euclidienne au lieu de matrice de covariance (ACP)

# Visualisation des données
# Heatmap
# e code ci-dessous permet de 1) calculer une distance entre chaque paire d’échantillon et de 2) tracer cette heatmap à partir de ces distances.
# t : transpose, dist : calcule distance, as.matrix : convertir en matrice
sampleDists <- as.matrix(dist(t(pseudoCounts)))
# scale = none : pas de normalisation, 
heatmap.2(sampleDists, scale = "none", col = bluered(100), trace = "none", density.info = "none", margins=c(6,6))

# Les couleurs chaudes comme le rouge indiquent des valeurs de distance plus faibles, ce qui suggère que les échantillons sont plus similaires en termes de profils d'expression des gènes.
# Les couleurs froides comme le bleu indiquent des valeurs de distance plus élevées, ce qui suggère que les échantillons sont plus différents ou divergents.

# On remarque que les échantillons basaux sont similaires en terme de profil d'expression de gènes
# A 0 : on a des profils identiques, raisonne en terme de distance (similarité)

# 3 questions à se poser : 
# Est ce que réplicat biologique similaire entre eux
# Est ce que agrégation cohérente avec 2 types de cellules
# Qu'est-ce qu'il passe qd on change condition souris

# Grosses diff entre basal et luminal.
# A, B : virgin
# C, D : pregnant
# E, F : lactate
# Entre pregnant et lactate moins grande distance qu'avec lactate et virgin
# Bcp plus changement et diff d'expression pr les luminales

plotMDS(pseudoCounts, col=as.numeric(SamplesInfo$condition))
# Duplicat super similaires entre eux super : séparation avant tout sur le type cellulaire, car axe 67% de la variabilité
# C'est conforme vis-à-vis des duplicats, pê que si on avait 3e dimension on verra mieux les distances (comme dans heatmap quoi)
# PAS PARLER DE DISTANCE LORS D'UNE ACP

# Identification des DEGs

dgeFull <- calcNormFactors(dgeFull, method="TMM")
dgeFull.group <- DGEList(CountData, remove.zeros = TRUE, group = dgeFull$samples$cellType)
dgeFull.group$samples$norm.factors <- dgeFull$samples$norm.factors
dgeFull.group <- estimateCommonDisp(dgeFull.group)
dgeFull.group <- estimateTagwiseDisp(dgeFull.group)

dgeExactTest <- exactTest(dgeFull.group)
resExactTest <- topTags(dgeExactTest, n = nrow(dgeExactTest$table))
ResultsDEGs.Full <- resExactTest$table

# À quoi correspondent les différentes colonnes de ce tableau ?

# logFC : Logarithme du rapport des expressions moyennes (fold change) entre les deux groupes comparés 
# (cellules luminales contre cellules basales). Une valeur positive indique une expression 
# plus élevée dans le premier groupe et une valeur négative indique une expression plus élevée dans 
# le second groupe (généralement le groupe "contrôle").

# logCPM : Logarithme du comptage moyen par million de reads. C'est une mesure de l'abondance générale du 
# gène dans l'ensemble des échantillons.

# PValue : La p-valeur de l'hypothèse nulle que le gène a le même niveau d'expression dans les deux groupes comparés. 
# Une p-valeur faible indique que les différences observées dans l'expression du gène entre les deux groupes sont
# statistiquement significatives.

# FDR : Le taux de fausse découverte (False Discovery Rate), ajusté pour les comparaisons multiples. Il 
# s'agit d'une correction appliquée aux p-valeurs pour tenir compte du problème des tests multiples, où 
# le test de plusieurs hypothèses simultanément augmente la probabilité d'obtenir des résultats faussement positifs.


# Pourquoi a-t-on une p-valeur par gène ? Pourquoi est-il important de faire une correction pour tests multiples ?

#Les p-valeurs par gène sont importantes car elles permettent de tester l'hypothèse nulle pour chaque gène 
# individuellement. Qd on effectue de nombreux tests (un pour chaque gène), la probabilité de trouver au moins 
# un résultat significatif par hasard augmente. La correction pour tests multiples (l'ajustement du FDR), 
# permet de contrôler cette probabilité et de réduire le nombre de découvertes faussement déclarées comme 
# significatives.

sum(ResultsDEGs.Full$FDR < 1e-4 & abs(resExactTest$table$logFC) > 3, na.rm=T) # log(2) = 1, diff au moins 2 fois supérieure
  
volcanoData <- cbind(resExactTest$table$logFC, -log10(resExactTest$table$FDR))
colnames(volcanoData) <- c("logFC", "-log10(p-value)")
DEGs <- resExactTest$table$FDR < 1e-4 & abs(resExactTest$table$logFC) > 3
point.col <- ifelse(DEGs, "red", "black")
plot(volcanoData, pch = 16, col = point.col, cex = 0.5)

# Différentiellement exprimés : gauche, bcp moins exprimés, à droite : bcp plus exprimés

# LogFc : différence des moyennes
# FDR : à quel pt mes ech sont dispersés (+ petit + confiance)
# Différence intensité et significativité
# S/ quels paramètres joués : gènes précis : on veut FDR petite
# Sp pas très bien séquencés : FDR plus importante, plus d'intra variabilité... c pg

# Croiser les informations entre NCBI et Ensemble
# Ensemble : regulatory build, info évolutive (quel autre gene orthologue)n gene expression etc.
# Sur NCBI : infos complt, 
# BDD UCSC : résultat classification fonctionnelle du génome

DiffExpressedGenes <- row.names(ResultsDEGs.Full[which(ResultsDEGs.Full$logFC>1 & ResultsDEGs.Full$FDR<0.05),])
write.table(file="Results_DEGs.txt", DiffExpressedGenes, quote=F, sep="\t", col.names=F, row.names=F)

