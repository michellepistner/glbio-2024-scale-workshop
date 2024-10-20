library(fgsea)
setwd("/home/hh/all_work/repos/glbio-2024-scale-workshop/Part_3_GSEA")

raw_lfc_data <- read.csv("thyroid_est_lfcs_songbird.tsv", sep="\t")
gmt.file <- gmtPathways("h.all.v7.4.symbols.gmt")

lfcs <- raw_lfc_data[,3]
names(lfcs) <- raw_lfc_data[,1]

fgsea.error <- function(lfcs, pathways, minSize=15, maxSize=500,
                        epsilons=c(-0.5,0,0.5)) {
    all_fgsea_res <- NULL
    for(epsilon in epsilons) {
        fgsea_res <- fgsea(stats=lfcs+epsilon, pathways=pathways,
                           minSize=minSize, maxSize=maxSize,
                           nPermSimple=5000)
        row.names(fgsea_res) <- NULL
        fgsea_res <- cbind(epsilon, fgsea_res)
        all_fgsea_res <- rbind(all_fgsea_res, fgsea_res)
    }
    return(all_fgsea_res)
}

names(gmt.file)

## Run vanilla fgsea
simple_fgsea_res <- fgsea(stats=lfcs, pathways=gmt.file, nPermSimple=5000)

## Run LFC Sensitivity Analysis FGSEA
lfc_fgsea_res <- fgsea.error(lfcs, gmt.file,
                             epsilon=c(-0.4, 0, 0.4))

simple_fgsea_res[(pathway=="HALLMARK_MYOGENESIS")|(pathway=="HALLMARK_INTERFERON_GAMMA_RESPONSE")]
lfc_fgsea_res[(pathway=="HALLMARK_KRAS_SIGNALING_DN")|(pathway=="HALLMARK_INFLAMMATORY_RESPONSE")]

