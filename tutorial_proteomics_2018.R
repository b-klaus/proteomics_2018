## ----options, include=FALSE----------------------------------------------
library(knitr)
options(digits = 3,
        width = 80)
opts_chunk$set(echo = TRUE, tidy = FALSE, include = TRUE,
               dev = 'png', fig.width = 6, fig.height = 3.5,
               comment = '  ', dpi = 300,
               cache = TRUE, warning = FALSE)

## ----required packages, echo = FALSE, warning=FALSE, results="hide"------
suppressPackageStartupMessages({
library(BiocStyle)
library(knitr)
library(Biobase)
library(geneplotter)
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(gplots)
library(RColorBrewer)
library(arrayQualityMetrics)
library(stringr)
library(matrixStats)
library(genefilter)
library(limma)
library(openxlsx)
library(smoothmest)
library(tidyr)
library(vsn)
library(MSnbase)
library(pheatmap)
library(fdrtool)
library(purrr)
library(tidyverse)
library(Hmisc)
library(gridExtra)
library(Peptides)
library("genefilter")
})

## ----required packages and data, echo = TRUE-----------------------------
library(BiocStyle)
library(knitr)
library(Biobase)
library(geneplotter)
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(gplots)
library(RColorBrewer)
library(arrayQualityMetrics)
library(stringr)
library(matrixStats)
library(genefilter)
library(limma)
library(openxlsx)
library(smoothmest)
library(tidyr)
library(vsn)
library(MSnbase)
library(pheatmap)
library(fdrtool)
library(purrr)
library(tidyverse)
library(Hmisc)
library(gridExtra)
library(Peptides)
library("genefilter")

glog2 <- function(x) ((asinh(x)-log(2))/log(2))

# fcuntion to retrieve annotation data from Uniprot

getUniprotGoodies  <- function(query) # columns)
{
    ## query and columns start as a character vectors
    qstring <- paste(query, collapse="+or+")
    #cstring <- paste(columns, collapse=",")
    uri <- 'http://www.uniprot.org/uniprot/?query='
    fullUri <- paste0(uri,qstring,'&format=tab')#&columns='),cstring)
    dat <- read.delim(fullUri, stringsAsFactors=FALSE)
    ## now remove things that were not in the specific original query...
    dat <- dat[dat[,1] %in% query,]
    dat
}


plotProtein <- function(Acc){
    
    tmp <- exprs(time_data_vsn)[Acc, ]
    pl <- (qplot(time, tmp, color = time, main = Acc) 
            + ylab("protein_exp_log2"))
    pl +  scale_color_brewer(type="qual", palette=2)
    }


# function to compute a robust mean 
robustMean <- function(x){
    if(length(x) == 1){return(x)}
    else{
        return(smhuber(x)$mu)
    } 
}


## ----impMeta-------------------------------------------------------------
expDesign <- as.data.frame(read_csv("metadata.csv"))
rownames(expDesign) <- expDesign$tmt.label

expDesign

## ----importData, eval = TRUE---------------------------------------------
load("peptide_data.RData")

peptide_data

peptide_data$fraction <- gsub("(.+R1_)|(.raw$)", "", peptide_data$source_file)
peptide_data$fraction <- factor(peptide_data$fraction,
                               ordered = TRUE,
                               levels = c("F01","F02","F03","F04","F05",
                                        "F06","F07","F08","F09","F10",
                                        "F11","F12"))

protein_data <- read_tsv("H0027_student_PEP_merged_results_20180202_1618_proteins.txt")

protein_data


## ----totalNoProts,  warning=FALSE----------------------------------------
ggplot(data = protein_data, aes(ifelse(qupm >= 2, "qupm >= 2", "qupm <2"))) +
  geom_bar() +
  xlab("") +
  geom_text(stat='count', aes(label = ..count..), nudge_y = 500) +
  ggtitle(paste("total number of identified proteins:", nrow(protein_data))) +
  coord_flip()

## ----protConc,  warning=FALSE--------------------------------------------
sub <- protein_data %>%
  dplyr::select(starts_with("signal_sum")) %>%
  gather()

ggplot(data = sub,aes(gsub("signal_sum_","",key), log2(value))) +
  geom_boxplot() +
  xlab("TMT-label") + ylab("log2(signal_sum)") +
  ggtitle("overview of log2(signal_sum) distributions")

rm(sub)

## ----pepMods-------------------------------------------------------------
mods <- peptide_data$modifications
mods <- gsub("[0-9]+;",";",mods)
mods <- gsub(";$","",mods)
mods <- unlist(strsplit(mods,split = "; "))
mods <- na.omit(mods)

qplot(mods,geom = "bar") +
  coord_flip() +
  ggtitle("total number of identified peptide modifications")

## ----trypticDigest-------------------------------------------------------
peptide_data <- peptide_data %>%
  mutate(sequence = toupper(sequence)) %>%
  mutate(no.of.miscleavages = gsub("[K,R]$","",sequence)) %>%
  mutate(no.of.miscleavages = gsub("[^KR]","",no.of.miscleavages)) %>%
  mutate(no.of.miscleavages = nchar(no.of.miscleavages))

bar_pl <- ggplot(data = peptide_data, aes(x = no.of.miscleavages)) +
  stat_count()

bar_data <- layer_data(bar_pl, 1) %>%
  mutate(prop_rounded = paste0(round(prop, 3)*100, "%"))

ggplot(data = bar_data, aes(x = x, y = y)) +
  geom_bar(stat = "identity") +
  geom_text(mapping = aes(x = x, y = y, label = prop_rounded), 
            nudge_y = 10000) +
  xlab("number of miscleavages") +
  ylab("number of peptides")


## ----peptideLengths------------------------------------------------------
ggplot(data = peptide_data, aes(nchar(sequence))) +
  geom_bar() +
  xlab("peptide lengths [aa]") +
  ylab("number of peptides") + 
  facet_wrap(~ifelse(in_quantification_of_protein == 1, 
                     "peptides used for protein quantification", 
                     "peptides not used for protein quantification"), 
             ncol = 1)

## ----TMTvsMS1, warning=FALSE---------------------------------------------
peptide_data$sig_sum <- peptide_data %>%
  dplyr::select(starts_with("sig")) %>%
  apply(1, sum, na.rm = TRUE)

scp <- ggplot(data = subset(peptide_data,in_quantification_of_protein == 1),
       mapping = aes(log2(peak_intensity), log2(sig_sum))) +
  stat_bin2d(binwidth = 0.2) +
  scale_fill_gradientn(colours = c("#377eb8","#984ea3","#e41a1c",
                                   "#ff7f00","#ffff33")) +
  xlab("log2(ms1-peak intensity)") +
  ylab("log2(sum of TMT-label signal)") 

scp + 
  annotate("text", label = paste0("Cor: ", 
                                  round(cor(layer_data(scp)$x,
                                            layer_data(scp)$y, 
                                            method = "spearman"), 2)),
           x = 22, y = 32, size = 5, colour = "coral3")

## ----plotMS1Fracs, fig.height=18, fig.width=15---------------------------
ggplot(data = subset(peptide_data, in_quantification_of_protein == 1),
       aes(`retention time`, peak_intensity)) +
  geom_line() +
  facet_wrap(~fraction)

## ----pepsPerFrac---------------------------------------------------------
ggplot(data = peptide_data, 
       aes(fraction, fill = ifelse(in_quantification_of_protein == 1,
                                   "quantified","not quantified"))) +
       geom_bar(position = position_dodge()) + 
       xlab("fraction") +
       ggtitle(paste("total number of identified peptides:", 
                     nrow(protein_data))) +
       scale_fill_brewer(palette = "Set1", name = "")

## ----C18column-----------------------------------------------------------
ggplot(data = subset(peptide_data,in_quantification_of_protein == 1),
       aes(`retention time`, fraction)) +
  stat_bin2d(binwidth = 1) +
  scale_fill_gradientn(colors = c("#377eb8","#984ea3","#e41a1c",
                                  "#ff7f00","#ffff33")) +
  xlab("retention time [min]")

## ----pIHeat--------------------------------------------------------------
peptide_data$pI <- Peptides::pI(peptide_data$sequence)
ggplot(data = subset(peptide_data,in_quantification_of_protein == 1),
       aes(`retention time`,fraction)) +
  stat_summary_2d(aes(z = pI),binwidth = 0.5) +
  scale_fill_gradientn(colours = c("#377eb8","#984ea3","#e41a1c",
                                   "#ff7f00","#ffff33"),
                       name = "isoelectric point pI") +
  xlab("retention time [min]") +
  ggtitle("prefractionation efficiency",
          subtitle = "by isoelectic point")

## ----hydrophobicityHeat--------------------------------------------------
peptide_data$gravy <- Peptides::hydrophobicity(peptide_data$sequence)
ggplot(data = subset(peptide_data, in_quantification_of_protein == 1),
       aes(`retention time`,fraction)) +
  stat_summary_2d(aes(z = gravy), binwidth = 0.5) +
  scale_fill_gradientn(colours = c("#377eb8","#984ea3","#e41a1c",
                                   "#ff7f00","#ffff33"),
                       name = "hydrophobicity (gravy score)") +
  xlab("retention time [min]") +
  ggtitle("prefractionation efficiency", 
          subtitle = "by hydrophobicity")

## ----fracGain------------------------------------------------------------
id_data <- NULL

pfractions <- levels(peptide_data$fraction)
for (i in seq_along(pfractions)) {
  sub <- subset(peptide_data,fraction %in% pfractions[1:i]&
                  in_quantification_of_protein==1)
  if (nrow(sub) > 0) {
    id_sub <- data.frame(fraction = pfractions[i],
                       unique.ids = length(unique(sub$protein_id)))
    id_data <- rbind(id_data,id_sub)
    rm(id_sub)
  }
}

id_data$fraction <- factor(id_data$fraction,
                         ordered = TRUE, levels = pfractions)
id_data$file <- "file"
ggplot(data = id_data, aes(fraction,unique.ids)) +
  geom_line(aes(group = file)) +
  geom_point() +
  ggtitle(label = "cum. sum of unique protein id's vs number of fractions")

## ----mascotScores--------------------------------------------------------
ggplot(data = peptide_data,aes(score,
       fill = factor(in_quantification_of_protein))) +
  geom_density(aes(y = ..count..), alpha = 0.5) +
  scale_fill_discrete(guide = guide_legend(title = "used for quant"))

## ----lengthVSscore-------------------------------------------------------
scp_mascot <- ggplot(data = subset(peptide_data, 
                   in_quantification_of_protein == 1), 
       aes(nchar(sequence),score)) +
  stat_bin2d(binwidth = 1) +
  geom_smooth(method = "lm") +
  scale_fill_gradientn(colours =
                         c("#377eb8","#984ea3",
                           "#e41a1c","#ff7f00","#ffff33")) +
  xlab("length of peptide [aa]")

scp_mascot + 
  annotate("text", label = paste0("Cor: ", 
                                  round(cor(layer_data(scp_mascot)$x,
                                            layer_data(scp_mascot)$y, 
                                            method = "spearman"), 2)),
           x = 15, y = 190, size = 5, colour = "coral3")

## ----uniquePepPerProt----------------------------------------------------
ggplot(data = subset(protein_data,qupm >= 2), aes(qupm)) +
  geom_bar()

## ----weightVSpeptides, warning=FALSE-------------------------------------
scp_mw <- ggplot(data = subset(protein_data, qupm >= 2),
       aes(as.numeric(as.character(mw)), upm)) +
  geom_point( alpha = I(0.1)) +
  geom_smooth(method = "loess") +
  scale_fill_gradientn(colours = c("#377eb8","#984ea3",
                                   "#e41a1c","#ff7f00","#ffff33")) +
  xlab("molecular weigth [Da]") +
  ylab("unique peptides per protein")
scp_mw + 
  annotate("text", label = paste0("Cor: ", 
                                  round(cor(layer_data(scp_mw)$x,
                                            layer_data(scp_mw)$y,
                                            use = "pairwise.complete.obs",
                                            method = "spearman"), 2)),
           x = 1e5, y = 190, size = 5, colour = "coral3")

## ----abunVSpeptides, warning=FALSE---------------------------------------
scp_top3 <- ggplot(data = subset(protein_data,qupm >= 2),aes(top3,upm)) +
  geom_point(alpha = I(0.1) ) +
  geom_smooth(method = "loess") +
  scale_fill_gradientn(colours = c("#377eb8","#984ea3","#e41a1c",
                                   "#ff7f00","#ffff33")) +
  xlab("top3 [abundance]") +
  ylab("unique peptides per protein")
scp_top3 + 
  annotate("text", label = paste0("Cor: ", 
                                  round(cor(layer_data(scp_top3)$x,
                                            layer_data(scp_top3)$y,
                                            use = "pairwise.complete.obs",
                                            method = "spearman"), 2)),
           x = 6.5, y = 190, size = 5, colour = "coral3")

## ----filterProt----------------------------------------------------------
protein_data <- protein_data %>%
  filter(!grepl("###",gene_name), qupm >= 2)

## ----getIntensityData----------------------------------------------------
prot_matrix <- dplyr::select(protein_data,
                             signal_sum_126:signal_sum_130H) %>%
               as.matrix()

colnames(prot_matrix) <- str_remove(colnames(prot_matrix), "signal_sum_")
rownames(prot_matrix) <- protein_data$protein_id 


stopifnot(colnames(prot_matrix)  == expDesign$tmt.label[-length(expDesign$tmt.label)])


 
to_exclude <- apply(prot_matrix, 1, function(x){any(!is.finite(x))})
table(to_exclude)

prot_matrix <- prot_matrix[!to_exclude, ]

head(prot_matrix)

## ----fData---------------------------------------------------------------
feature_anno <- as.matrix(dplyr::select(protein_data, protein_id, gene_name,
                              top3, qupm, description))

rownames(feature_anno) <- protein_data$protein_id

feature_anno <- as.data.frame(feature_anno)[!to_exclude, ]

## ----sumexp, echo=FALSE, fig.show="asis"---------------------------------
par(mar=c(0,0,0,0))
plot(1,1,xlim=c(0,100),ylim=c(0,100),bty="n",
     type="n",xlab="",ylab="",xaxt="n",yaxt="n")
polygon(c(45,80,80,45),c(10,10,70,70),col=rgb(1,0,0,.5),border=NA)
polygon(c(45,80,80,45),c(68,68,70,70),col=rgb(1,0,0,.5),border=NA)
text(62.5,40,"assay(s)", cex = 1)
text(62.5,30,"e.g. 'exprs'", cex = 1)
polygon(c(20,40,40,20),c(10,10,70,70),col=rgb(0,0,1,.5),border=NA)
polygon(c(20,40,40,20),c(68,68,70,70),col=rgb(0,0,1,.5),border=NA)
text(30,40,"featureData", cex = 1)
polygon(c(45,80,80,45),c(75,75,90,90),col=rgb(.5,0,.5,.5),border=NA)
polygon(c(45,47,47,45),c(75,75,90,90),col=rgb(.5,0,.5,.5),border=NA)
text(62.5,82.5,"phenoData", cex = 1)

## ----createEset, dependson="summarizeData", eval = TRUE------------------

eset_raw <- ExpressionSet(prot_matrix, 
                   phenoData = AnnotatedDataFrame(expDesign[-nrow(expDesign),]),
                   featureData = AnnotatedDataFrame(feature_anno))
  
validObject(eset_raw)
pData(eset_raw)

# prots with very low expression / zero expression 
strange_prots <- apply(exprs(eset_raw), 1, function(x){min(x) < 1})
eset_raw <- eset_raw[!strange_prots, ]
save(eset_raw, file = "eset_raw.RData")


## ----arrayQualityMetricsRaw, eval = FALSE--------------------------------
## try(arrayQualityMetrics(expressionset = eset_raw,
##       outdir = "Report_for_eset_raw",
##     force = TRUE,
##     do.logtransform = TRUE,
##    intgroup = c("condition", "rep")
## ))
## 

## ----boxplotRaw----------------------------------------------------------
oligo::boxplot(eset_raw, transfo = log2, range = 0)

## ----vsnNorm, dependson="createEset"-------------------------------------
vsn_fit <- vsn2(eset_raw)
eset_vsn <- predict(vsn_fit, newdata = eset_raw)

oligo::boxplot(eset_vsn, transfo = identity)

save(eset_vsn, file = "eset_vsn.RData")

## ----arrayQualityMetricsAfterVSN, eval = FALSE---------------------------
## try(arrayQualityMetrics(expressionset = eset_vsn,
##       outdir = "Report_for_eset_after_vsn",
##     force = TRUE,
##     do.logtransform = FALSE,
##    intgroup = c("condition", "rep")
## ))
## 

## ----topVarPCA, dependson="vsn_norm", eval = TRUE------------------------

ntop <- 500

Pvars <- rowVars(exprs(eset_vsn))
select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, 
        length(Pvars)))]
PCA <- prcomp(t(exprs(eset_vsn)[select, ]), scale = TRUE)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)


dataGG = data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], 
                    PC3 = PCA$x[,3], PC4 = PCA$x[,4], 
                    replicate = pData(eset_vsn)$rep, 
                    condition = pData(eset_vsn)$condition)
                  
(qplot(PC1, PC2, data = dataGG, color =  condition, 
       main = "PC1 vs PC2, top variable proteins", size = I(6), 
       shape = replicate)
 + labs(x = paste0("PC1, VarExp:", round(percentVar[1],4)),
        y = paste0("PC2, VarExp:", round(percentVar[2],4)))
 + scale_colour_brewer(type = "qual", palette = 2)
 )



## ----createDesign, eval = TRUE-------------------------------------------
condition <- as_factor(pData(eset_vsn)$condition )
des <- model.matrix(~ condition)
colnames(des) <- str_remove(colnames(des), "condition")
des

## ----compare_to_time0, dependson="createDesign", eval = TRUE-------------

lm_fit <- eBayes(lmFit(eset_vsn, design = des))
            
limma_table <- topTable(lm_fit, sort.by = "t",  coef = "single_ko",
                        number = Inf)

as.tibble(limma_table)


table(limma_table$adj.P.Val < 0.0001)

## ----pvalHist------------------------------------------------------------
hist(limma_table$P.Value, col = "lavender")

## ----volcPlot------------------------------------------------------------
ggplot(data = limma_table,aes(logFC,-log10(P.Value))) +
  geom_point(alpha = I(0.5)) +
  geom_vline(xintercept = 0) +
  xlab("log2(fold change)") +
  geom_text(data = subset(limma_table, adj.P.Val <= 0.0001), 
            aes(label = gene_name),
            size = 2,nudge_y = 0.15,check_overlap = TRUE) +
  ggtitle("volcano plot", subtitle = "single_ko vs wt")

## ----MAplot--------------------------------------------------------------
ggplot(data = limma_table,aes(AveExpr, logFC)) +
  geom_point(alpha = I(0.2)) +
  geom_hline(yintercept = 0) +
  ylab("log2(fold change)") +
  geom_text(data = subset(limma_table,adj.P.Val <= 0.0001), 
            aes(label = gene_name),
            size = 2, nudge_y = 0.08,check_overlap = TRUE) +
  ggtitle("MA plot",subtitle = "single_ko vs wt")

## ----top3VSlogFC---------------------------------------------------------
ggplot(data = limma_table, 
       aes(logFC,as.numeric(as.character(top3)))) +
  geom_point(alpha = I(0.2)) +
  geom_vline(xintercept = 0) +
  xlab("log2(fold change)") +
  ylab("top3 [abundance]") +
  geom_text(data = subset(limma_table,adj.P.Val <= 0.0001),
            aes(label = gene_name),
            size = 2, nudge_y = 0.08,check_overlap = TRUE) +
  ggtitle("top3 plot",subtitle = "single_ko vs wt")

## ----lossOfFDRControl----------------------------------------------------
n <- 20
p <- 10000

true_fcs <- numeric(10000)
idx_ab <- sample(x = p, 1000)
true_fcs[idx_ab] <- rnorm(length(idx_ab), mean = 2, sd = 0.5)

x <- matrix(rnorm(n*p), nrow = p, ncol = n)
fac <- factor(c(rep(0, 10), rep(1, 10)))
x[idx_ab, fac == 1] <- x[idx_ab, fac == 1]  + true_fcs[idx_ab]
rt <- rowttests(x, fac)

## ----ordinaryFDR---------------------------------------------------------
FDR <- p.adjust(rt$p.value, method = "BH")

P <- which(FDR < 0.1)
TP <- intersect(idx_ab, P)
FP <- setdiff(P, idx_ab)

FDP <- length(FP) / length(P)

FDP

## ----postHocFDR----------------------------------------------------------
fcs <- rt$dm

thresh_hits <- which(abs(fcs) > 2 & FDR < 0.1)

P_thresh <- thresh_hits 

true_fcs_2 <- intersect(idx_ab,  which(abs(true_fcs) > 2))

FP_thresh <- setdiff(P_thresh, true_fcs_2)

FDP_thresh <- length(FP_thresh) / length(P_thresh)

FDP_thresh

## ----treatFit------------------------------------------------------------

treat_res <- topTreat(treat(lmFit(x, design = model.matrix(~ fac)), lfc = 2),
                       coef = 2, number = Inf, sort.by = "none")

P_treat <- which(treat_res$adj.P.Val < 0.1)

FP_treat <- setdiff(P_treat, true_fcs_2)

FDP_treat <-  length(FP_treat) / length(P_treat)

FDP_treat

## ----seesionInfo, results='markup'---------------------------------------
sessionInfo()

