# Pa_biofilm_BIOLOG

## makeTidyTables

A script to analyze crystal violet data obtained from experiments on <i>Pseudomonas aeruginosa</i>. The script takes a table of data containing information on the OD (optical density) and CV (crystal violet) measurements of <i>P. aeruginosa</i> grown under different nutrient conditions and oxygen availability, then plots the joint distribution of CV and OD measurements for anaerobic and aerobic conditions, computes CV and OD values corrected for negative control using fitlm, Finally it uses fitlm to determine the change in biofilm formation when switching from aerobic to anaerobic conditions for each nutrient. The final plot shows the results of this analysis for nutrients where there is evidence they can be utilized for growth above negative control, either aerobically or anaerobically.
