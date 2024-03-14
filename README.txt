Replication files for the paper "Identification and Inference of Network Formation Games with Misclassified Links" by Luis E. Candelaria and Takuya Ura.

The master file is named 'main_MLE.m'. Prior to running the code, 
1. Change the [PATH] in the file 'main_MLE.m' to the actual directory path of this folder.
2. Create the following folders in the current directory:
'code', 'raw\_data', 'directed\_adjacency\_matrices', 'dMU', 'lendmoney\_graphs', and 'results\_mis'.

The structure of this directory is the following:

1. 'raw\_data' folder:
Download data from
http://www.stanford.edu/~jacksonm/IndianVillagesDataFiles.zip. 
Unzip and put the contents of the folder '2010-0760\_Data' directly into this folder.

2. 'code' folder:
Download and insert the Matlab functions "matrix2latex" and "Multiprod" into this folder. Functions are available here 
https://uk.mathworks.com/matlabcentral/fileexchange/4894-matrix2latex
https://uk.mathworks.com/matlabcentral/fileexchange/8773-multiple-matrix-multiplications-with-array-expansion-enabled

3. 'directed\_adjacency\_matrices' folder: 
Contains csvs files of 'lendmoney' and 'rel' generated from 
'D100_drop_notype_all.m'and 'D100_gendirected.m.'

4. 'dMU' folder:
Contains dyad-level covariates and nonparametric first-stage estimates generated from 'D200_dMU.m' and 'G100_PsiGamma.m'.

5. 'results_mis' folder:
Contains the results from the 'main_MLE.m' and 'V100_Avarmis.m'

Acknowledgements:

We benefited substantially from the data and code used in the following papers.

* Leung, M. P. (2015). Two-step estimation of network-formation models with incomplete information. Journal of Econometrics, 188(1), 182-195. 
* Abhijit Banerjee, Arun Chandrasekhar, Esther Duflo, Matthew O. Jackson (2013) ``The Diffusion of Microfinance,'' _Science_, 341 (6144)
* Matthew O. Jackson, Tomas Rodriguez-Barraquer, Xu Tan (2012) ``Social Capital and Social Quilts: Network Patterns of Favor Exchange,'' _American Economic Review_, 102 (5).
