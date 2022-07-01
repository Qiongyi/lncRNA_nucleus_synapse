
####################################################################################################################
############                        1. lncRNA Capture-Seq data analysis                                 ############
####################################################################################################################

### 1) Use cutadapt (v1.17) to clip low-quality nucleotides and adaptor sequences
# Bases lower than a defined Phred quality threshold (default: 20) at the 3′ end were trimmed off from each read using cutadapt (http://code.google.com/p/cutadapt/). Next, known Illumina primers and adaptor sequences were clipped off from each read by cutadapt, which computes sensitive semi-global alignments of all the reads against all the primer/adaptor sequences, allowing gapped and mismatched alignments.

pipeline_cutadapt.pl sample_info.txt /clusterdata/uqqzhao/illumina/202008_WeiSiang_CapSeq_BredyLab/ori /clusterdata/uqqzhao/illumina/202008_WeiSiang_CapSeq_BredyLab/lncRNACapSeq/cutadapt


### 2) Processed reads were aligned against the ribosome RNA and PhiX reference sequences using bowtie2 (version 2.3.4.2).
indir=/clusterdata/uqqzhao/illumina/202008_WeiSiang_CapSeq_BredyLab/lncRNACapSeq
pipeline_alignment.pl bowtie2 $indir/config.txt /illumina/reference/mm10/mm10_rRNA $indir/cutadapt $indir/rRNA_bowtie2
pipeline_alignment.pl bowtie2 $indir/config.txt /illumina/reference/PhiX/PhiX $indir/cutadapt $indir/PhiX_bowtie2



### 3) mapping against the mouse genome (mm10) using HISAT2 (v2.1.0); convert “SAM” files to “BAM” files, remove duplicate reads, sort and index the “BAM” files. 
# To avoid the artefact signals potentially introduced by misalignments, we only kept properly PE aligned reads with mapping quality at least 20 for downstream analyses.

cd $indir
pipeline_alignment.pl hisat2 config.txt /illumina/reference/mm10/HISAT2_index/mm10 ./cutadapt ./HISAT2


### 4) run StringTie - three rounds of StringTie (v2.1.4) was applied to i) perform reference-guided transcriptome assembly by supplying the GENCODE annotation file (V25) with the “-G” option for each sample, ii) generate a non-redundant set of transcripts using the StringTie merge mode, and iii) quantifying the transcript-level expression for each sample, with the option of “-e -G merged.gtf”. 

# 4.1) perform reference-guided transcriptome assembly for each sample
for i in "synapse1" "synapse2" "synapse3" "synapse4"
do
/clusterdata/uqqzhao/tools/StringTie/stringtie-2.1.4.Linux_x86_64/stringtie ./HISAT2/$i.rmdup.Q20.sort.bam -p 2 -o $i.gtf -e -b ./StringTie/$i.ctab -G /illumina/reference/GENCODE/gencode.vM25.annotation.gtf
done

# 4.2) generate a non-redundant set of transcripts using the StringTie merge mode. Note we also merged with previous nucleus samples so that we can directly perform the comparison analyses between these two compartments.

nuc_dir=/illumina/Data/others/201912_Wei_Nuclear_LncRNA_CapSeq_BredyLab/StringTie/StringTie_2.1.4

/clusterdata/uqqzhao/tools/StringTie/stringtie-2.1.4.Linux_x86_64/stringtie --merge -G /clusterdata/uqqzhao/reference/GENCODE/gencode.vM25.annotation.gtf -o merge_with_genecode.gtf synapse1.gtf synapse2.gtf synapse3.gtf synapse4.gtf $nuc_dir/EXT1.gtf $nuc_dir/EXT2.gtf $nuc_dir/EXT3.gtf $nuc_dir/RC4.gtf $nuc_dir/RC5.gtf $nuc_dir/RC6.gtf


# 4.3) quantitate the transcript-level expression for each sample, with the option of “-e -G merged.gtf”
for i in "synapse1" "synapse2" "synapse3" "synapse4"
do
/clusterdata/uqqzhao/tools/StringTie/stringtie-2.1.4.Linux_x86_64/stringtie ./HISAT2/$i.rmdup.Q20.sort.bam -p 2 -o $i.eG.gtf -e -b $i.ctab -G merge.gtf &
done

nuc_dir=/illumina/Data/others/201912_Wei_Nuclear_LncRNA_CapSeq_BredyLab/HISAT2
for i in "EXT1" "EXT2" "EXT3" "RC1" "RC2" "RC3"
do
/clusterdata/uqqzhao/tools/StringTie/stringtie-2.1.4.Linux_x86_64/stringtie $nuc_dir/$i.rmdup.Q20.sort.bam -p 2 -o $i.eG.gtf -e -b $i.ctab -G merge.gtf &
done

# clean StringTie GTF (remove transcripts that are not listed in "merge.gtf" - note this is a potential bug that sometimes exists in the StringTie output file.)
for i in "synapse1" "synapse2" "synapse3" "synapse4" "EXT1" "EXT2" "EXT3" "RC1" "RC2" "RC3"
do
StringTie_clean_GTF.pl merge.gtf $i.eG.gtf $i.clean.gtf
done


# extract FPKM for all samples
ls -l EXT?.clean.gtf >group1
ls -l RC?.clean.gtf >>group1
ls -l syn*.clean.gtf >group2
StringTie_GTF2FPKM.pl group1 group2 nucleus_synapse_FPKM.xls



### 5) run ballgown - see ballgown.R


### 6) Extract nucleus enriched lncRNA list and synapse enriched lncRNA list
# 6.1) remove known protein-coding transcripts based on GENCODE annotation
Remove_known_coders.pl /clusterdata/uqqzhao/reference/GENCODE/gencode.vM25.annotation.gtf Nucleus.vs.Synapse.lncRNACapture.ballgown.xls others.xls protein_coding.xls

# 6.2) filter lowly expressed transcripts, FPKM>1 in at least one group
awk '$1=="gene_id" || ($NF>1 || $(NF-1)>1)' nucleus_synapse_FPKM.xls > nucleus_synapse_FPKM.filter.xls

# 6.3) apply both 6.1 & 6.2 (remove protein-coding, remove lowly expressed transcripts)
Extract_known_lnRNAs.pl others.xls nucleus_synapse_FPKM.filter.xls nucleus_synapse_FPKM.qval.known.xls nucleus_synapse_FPKM.qval.xls
awk '$1=="gene_id" || ($NF<0.05 && $(NF-2)>2 && $(NF-4)>$(NF-3))' nucleus_synapse_FPKM.qval.xls >  lncRNA_nucleus_enriched.xls
awk '$1=="gene_id" || ($NF<0.05 && $(NF-2)<0.5 && $(NF-4)<$(NF-3))' nucleus_synapse_FPKM.qval.xls > lncRNA_synapse_enriched.xls


### 7) alternative splicing analysis using SUPPA (https://github.com/comprna/SUPPA)
python /scratch/qbi/uqqzhao/tools/SUPPA/SUPPA-2.3/suppa.py generateEvents -i merge_with_genecode.gtf -o merge -f ioi 
python /scratch/qbi/uqqzhao/tools/SUPPA/SUPPA-2.3/suppa.py generateEvents -i merge_with_genecode.gtf -o merge -f ioe -e SE SS MX RI FL

lncRNA_AS_stats.pl lncRNA_nucleus_enriched.SINE.Len.xls merge lncRNA_nucleus_enriched.AS.stats.xls lncRNA_nucleus_enriched.AS.xls 
lncRNA_AS_stats.pl lncRNA_synapse_enriched.SINE.Len.xls merge lncRNA_synapse_enriched.AS.stats.xls lncRNA_synapse_enriched.AS.xls
