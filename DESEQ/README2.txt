RUN THE Rscript HMAP&PCA.R like this

./HMAP&PCA.R 

Dependencies: ALl pairwise .csv files, DESEQ package installed.

OUTPUT:Generates the PCA plots, Heatmaps for pairwise samples 



Alternatively, one may choose to run a similar script HMAP_PCA.R as:

./HMAP_PCA.R <input_file.csv>

Incase <input_file.csv> is missing, the program would print an error message and exit
For e.g.: ./HMAP_PCA.R

Error: At least one argument must be supplied (input file).n
Execution halted

output:Gnerates the PCA plots , Heatmaps for ALL VALID PAIRWISE SAMPLES provided as an argument
