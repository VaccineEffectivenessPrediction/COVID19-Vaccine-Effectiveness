# COVID19-Vaccine-Effectiveness

R codes for the article "Rapid evaluation of COVID-19 vaccine effectiveness against SARS-CoV-2 variants by genetic distance"

This repository includes the following three folders:
  Step 1 Preparing data;
  Step 2 Calculating Genetic Distance;
  Step 3 Data Summary and Create Figures.
Data analysis involved in this work can be performed with these codes. 

Step 1 Preparing data: 
The R file "PrepareProteinSequences.R" can be used to clean sequencing data and export data as csv fromat.  

Step 2 Calculating Genetic Distance: 
This folder includes two R files. The R file "CalculateGeneticDistanceofSprotein.R" can be used to calculate genetic distance on NTD, RBD and entire Spike protein. The R file "CalculateGeneticDistanceofnonSprotein.R" can be used to calculate genetic distance on non-S proteins of SARS-CoV-2 genome.

Step 3 Data Summary and Create Figures:
This folder includes five R files to create figures of the main text and extended data. The R file "Figure_1andExtendedDataFigure_1.R" is needed to compare distributions of vaccine efficacy or effectiveness (VE) and genetic distance for NTD, RBD and Spike protein. The R file "Figure_2andExtendedDataFigures_2-4.R" can be used to make scatter plot for different platforms and different vaccine products based on mismatch on a given genomic region to investigate the association between VE and genetic distance. The R file "Figure_3.R" is the code for the training-validation setting and "Figure_4.R" for predicting VE against variants and for California. The R file "Figure_5andExtendedDataFigure_5.R" is the script to generate a heatmap for exploring the possibility of region-specific vaccines. 

The PDF file "Supplementary Acknowledgement Table.pdf" reports accussion numbers of SARS-CoV-2 sequences. Viral sequence data were downloaded from the global initiative on sharing all influenza data (GISAID) at http://platform.gisaid.org/. We thank the contributions of all the health care workers and scientists, the GISAID team, and the submitting and the originating laboratories. 
