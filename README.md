# Pa_biofilm_BIOLOG

## makeTidyTables

A script to analyze crystal violet data obtained from experiments on <i>Pseudomonas aeruginosa</i>. The script takes a table of data containing information on the OD (optical density) and CV (crystal violet) measurements of <i>P. aeruginosa</i> grown under different nutrient conditions and oxygen availability. It uses fitlme to determine the change in biofilm formation when switching from aerobic to anaerobic conditions for each nutrient. Saves the tables in folder 'tidy_tables'.

## unsupervisedAnalysis

The script loads a tidy table of CV and OD data, plots the distributions of the data, and then generates a clustergram and a PCoA plot based on correlation distance. The clustergram can be customized to plot only the nutrients with evidence of growth, or all 190 nutrients. The PCoA plot shows the log2-fold change in biofilm formation in anaerobic vs. aerobic conditions for each nutrient.
