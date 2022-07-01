
### ballgown
setwd("/clusterdata/uqqzhao/illumina/202008_WeiSiang_CapSeq_BredyLab/lncRNACapSeq/StringTie_GencodeV25")

# on inode2
# module load R/4.0.3
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ballgown")

library(ballgown)

library(dplyr) # for the arrange function

#bg = ballgown(samples=c("RC4.ctab", "RC5.ctab", "RC6.ctab", "EXT1.ctab", "EXT2.ctab", "EXT3.ctab", "SynRC1.ctab", "SynRC3.ctab", "SynEXT1.ctab"), meas='all')
bg = ballgown(samples=c("RC4.ctab", "RC5.ctab", "RC6.ctab", "EXT1.ctab", "EXT2.ctab", "EXT3.ctab", "synapse1.ctab", "synapse2.ctab", "synapse3.ctab", "synapse4.ctab"), meas='all')
pData(bg) = data.frame(id=sampleNames(bg), group=c(1,1,1,1,1,1,0,0,0,0))

stat_results = stattest(bg, feature='transcript', meas='FPKM', getFC=TRUE, covariate='group')

#results_transcripts = data.frame(transcriptIDs=ballgown::transcriptNames(bg), geneIDs=ballgown::geneIDs(bg), geneNames=ballgown::geneNames(bg), stat_results)
results_transcripts = data.frame(transcriptIDs=ballgown::transcriptNames(bg), geneIDs=ballgown::geneIDs(bg), geneNames=ballgown::geneNames(bg), stat_results[,3:5])

results_transcripts = arrange(results_transcripts,pval)
head(results_transcripts)

write.table(results_transcripts,"Nucleus.vs.Synapse.lncRNACapture.ballgown.xls", row.names=FALSE, quote = FALSE, sep = "\t")

# save the sessionInfo
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")


