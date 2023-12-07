# TCR-repertoire-Machine-learning
TCR repertoire characteristics predict clinical response to adoptive CTL therapy against nasopharyngeal carcinoma 
Interested in reading paper, here is the link - https://pubmed.ncbi.nlm.nih.gov/34377592/

#=== Input and output files are not made available.

EBV TCR analysis

Step1
	Aim - Read Mixcr TCR clones input files and the associated phenotype file
	Script - Rscript_1_ReadInFiles.R
	Input  - /data/*_clones_full.txt,
		 phenotype_hla.txt
	
Step2
	Aim - Calculate Jaccard Amino acid similarity between samples and count the HLA shared between samples.
	Script - Rscript_2_Similarity.R
	Input - Step1
	Output - similarity.hla.RData

Step3
	Aim - Calculate Shannon entropy (diversity) values for Full TCR and CDR3 for each sample.
	Script - Rscript_3_Diversity.R
	Input - Step1
	Output - diversity.hla.RData
	
Step4
	Aim - Generate Figures 2A to 2F
	Script - Rscript_4_Figure_2Ato2F.R
	Input - similarity.hla.RData
	Output - 2A.survival.similarity.byResp.full.tcr.svg,
		 2B.survival.similarity.byResp.cdr3.svg,
		 2C.survival.similarity.byResp.full.tcr.svg,
		 2D.perclones.similarity.byResp.cdr3.svg,
		 2E.survival.hlaShared.byResp.svg,
		 2F.hla.shared.similarity.byResp.all.clones.cdr3.svg

Step5
	Aim - Generate Figures 1A to 1C
	Script - Rscript_5_Figure_1Ato1C.R
	Input - diversity.hla.RData
	Output - 1A.diversity.cdr3.svg,
		 1B.diversity.full.tcr.svg,
		 1C.diversity.hlaA.alleles.full.tcr.svg

Step6
	Aim - Generate Supplementary Figures 1A and 1B
	Script - Rscript_6_SuppFigure_1A_1B.R
	Input - diversity.hla.RData
	Output - Supp.1A.survival.diversity.cdr3.svg,
		 Supp.1B.survival.diversity.full.tcr.svg

Step7
	Aim - Generate Supplementary Figures 2A to 2E
	Script - Rscript_7_SuppFigure_2Ato2E.R
	Input - diversity.hla.RData
	Output - Supp.2A.diversity.hlaA.full.tcr.svg,
		 Supp.2B.diversity.hlaB.alleles.full.tcr.svg,
		 Supp.2C.diversity.hlaB.full.tcr.svg,
		 Supp.2D.diversity.hlaC.alleles.full.tcr.svg,
		 Supp.2E.diversity.hlaC.full.tcr.svg
		 

EBV MOTIF ANALYSIS


Step8
	Aim - Read Mixcr TCR clones input files and the associated phenotype file. Apply important filters to QC the clones.
	Script - Rscript_8_ReadInFiles_QC.R
	Input  - /data/*_clones_full.txt,
		 phenotype_hla.txt
	
Step9
	Aim - Calculate Motif cdr3 counts and accumulative cdr3 counts. Calculate Correlation between Response samples and Overall survival.
	Script - Rscript_9_MotifCounts_Correlation.R
	Input - Step8
	Output - RtoNRratio_Rcorrelation.RData,
		  motif_cdr3_counts.RData,
		  motif_accu_cdr3_counts.RData

Step10
	Aim - Run Machine learning Random forest to find best cutoff to select high accuracy motifs. Figure 3A.
	Script - Rscript_10_Figure_3A.R
	Input - RtoNRratio_Rcorrelation.RData,
		motif_cdr3_counts.RData
	Output - 3A.Accuracy.vs.RtoNR.cdr3.cutoff.svg
	
Step11
	Aim - R/NR Ratio vs Overall survival plot, unsupervised clustering, final motifs accuracy and Gini plots. Figures 3B 3C and 3D.
	Script - Rscript_11_Figure_3Bto3D.R
	Input - RtoNRratio_Rcorrelation.RData,
		motif_cdr3_counts.RData
	Output - 3B.Correlation.vs.RtoNR.cdr3.ratio.svg,
		 3C.dendrogram.svg,
		 3D.uniqueCDR3.meanDecreaseAccuracy.svg,
		 3D.uniqueCDR3.meanDecreaseGini.svg
		  

Step12
	Aim - ROC curve for SPDQ and SPFS. Figure 3F.
	Script - Rscript_12_Figure_3F.R
	Input - motif_cdr3_counts.RData
	Output - 3F.SPDQ.ROC.unique.cdr3.counts.svg,
		 3F.SPFS.ROC.unique.cdr3.counts.svg

Step13
	Aim - Boxplot (with response) and scatterplot (with overall survival) for SPDQ and SPFS. Figures 3E 3G 3H.
	Script - Rscript_13_Figure_3E_3G_3H.R
	Input - motif_accu_cdr3_counts.RData
	Output - 3E.SPDQ.resp.vs.logclonecount.svg,
		 3G.SPDQ.OS.vs.logclonecount.svg,
		 3H.SPDQ.hla.alleles.vs.logclonecount.svg,
		 3E.SPFS.resp.vs.logclonecount.svg,
		 3G.SPFS.OS.vs.logclonecount.svg,
		 3H.SPFS.hla.alleles.vs.logclonecount.svg

Step14
	Aim - Boxplot (with response) and scatterplot (with overall survival) for remaining motifs. Supplementary Figures 3A and 3B
	Script - Rscript_14_Supp_Figure_3A_3B.R
	Input - motif_accu_cdr3_counts.RData
	Output - Supp.3A.MOTIF.resp.vs.logclonecount.svg,
		 Supp.3B.MOTIF.OS.vs.logclonecount.svg
		 
