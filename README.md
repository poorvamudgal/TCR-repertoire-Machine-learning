# TCR-repertoire-Machine-learning
TCR repertoire characteristics predict clinical response to adoptive CTL therapy against nasopharyngeal carcinoma.

Interested in reading paper, here is the link - https://pubmed.ncbi.nlm.nih.gov/34377592/

#=== Input and output files are not made available.

EBV TCR analysis

Step1
	Aim - Read Mixcr TCR clones input files and the associated phenotype file -> 
	Script - Rscript_1_ReadInFiles.R
	
Step2
	Aim - Calculate Jaccard Amino acid similarity between samples and count the HLA shared between samples -> 
	Script - Rscript_2_Similarity.R

Step3
	Aim - Calculate Shannon entropy (diversity) values for Full TCR and CDR3 for each sample -> 
	Script - Rscript_3_Diversity.R
	
Step4
	Aim - Generate Figures 2A to 2F -> 
	Script - Rscript_4_Figure_2Ato2F.R
 
Step5
	Aim - Generate Figures 1A to 1C -> 
	Script - Rscript_5_Figure_1Ato1C.R
 
Step6
	Aim - Generate Supplementary Figures 1A and 1B -> 
	Script - Rscript_6_SuppFigure_1A_1B.R

Step7
	Aim - Generate Supplementary Figures 2A to 2E -> 
	Script - Rscript_7_SuppFigure_2Ato2E.R

EBV MOTIF ANALYSIS


Step8
	Aim - Read Mixcr TCR clones input files and the associated phenotype file. Apply important filters to QC the clones -> 
	Script - Rscript_8_ReadInFiles_QC.R

Step9
	Aim - Calculate Motif cdr3 counts and accumulative cdr3 counts. Calculate Correlation between Response samples and Overall survival -> 
	Script - Rscript_9_MotifCounts_Correlation.R

Step10
	Aim - Run Machine learning Random forest to find best cutoff to select high accuracy motifs. Figure 3A -> 
	Script - Rscript_10_Figure_3A.R

Step11
	Aim - R/NR Ratio vs Overall survival plot, unsupervised clustering, final motifs accuracy and Gini plots. Figures 3B 3C and 3D -> 
	Script - Rscript_11_Figure_3Bto3D.R

Step12
	Aim - ROC curve for SPDQ and SPFS. Figure 3F -> 
	Script - Rscript_12_Figure_3F.R

Step13
	Aim - Boxplot (with response) and scatterplot (with overall survival) for SPDQ and SPFS. Figures 3E 3G 3H -> 
	Script - Rscript_13_Figure_3E_3G_3H.R

Step14
	Aim - Boxplot (with response) and scatterplot (with overall survival) for remaining motifs. Supplementary Figures 3A and 3B -> 
	Script - Rscript_14_Supp_Figure_3A_3B.R
	 
