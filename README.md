# SpacePhenotyper
SpacePhenotyper: A computational method for translating information between different (but related) types of data
#### Pre-required installations before running SpacePhenotyper
Python libraries pandas, numpy, seaborn, scipy, matplotlib and sklearn are prerequired to be installed before running SpacePhenotyper
#### Input Data to SpacePhenotyper
Quantity vector of a cancer-related phenotype across patents, corresponding bulk gene expression data, and spatially resolved transcriptomic (SRT) data as cvs files
```
BulkExpression = pd.read_csv('Bulk-Expr-BRCA-TCGA.csv', index_col=0)           # Load the bulk gene expression matrix (columns are genes and rows are patients)
PhenotypeVector = pd.read_csv('PhenotypeVector-TP-BRCA-TCGA.csv', index_col=0) # Load the vector of phenotype quantity matched with the rows of bulk gene expression matrix (a vector across patients)
DataSRT = pd.read_csv('SRT-Expr-Breast-1.1.0-Smoothed.csv', index_col=0)       # SRT gene expression matrix (rows are genes, and columns are spots)
MetaData = pd.read_csv('SRT-tissue-positions-Breast-1.1.0.csv', index_col=0)   # Load spot locations in SRT data
```
#### To Run SpacePhenotyper
SpacePhenotyper takes a phenotype quantity vector and a bulk gene expression matrix across multiple patients, and SRT gene expression matrix and spot locations matched ascross spots  
```
 [EigenGene, EigenPatient, Result]= SpacePhenotyper(PhenotypeVector, BulkExpression, DataSRT, MetaData)  # Predict the phenotype quantity on spots of tissue slice  
```
#### Output of SpacePhenotyper
As result, returns the Eigen-Gene vector, Eigen-Patient vector, Cosine similarity vector (Predicted phenotype quantity on spots) and Scatter plot for the predicted phenotype quantity on the spot locations  
```
EigenGene.to_csv('Eigen-Gene.csv')                  # Save the Eigen-Gene vector as a cvs file
EigenPatient.to_csv('Eigen-Patient.csv')            # Save the Eigen-Patient vector as a cvs file
Result['SpacePhenotyper'].to_csv('PredictedValues-On-Spots.csv')  # Save the predicted phenotype quantity vector as a cvs file
plt.savefig(Method+'-Prediction-Plot.png', dpi = 300, bbox_inches = 'tight')    # Save the prediction plot for phenotype quantity
```
#### Python [Package](code)
* The SpacePhenotyper method is implemented in python and the codes are available as [Python Code](code/SpacePhenotyper.py) and [Jupyter Notebook](code/SpacePhenotyper.ipynb) modules.

#### Bulk data sets [Data Bulk](Bulk-data)
The phenotype quantity of patients (Tumor-Purity, Hazard, RCB, Immune and Proliferation) and bulk gene expression matrices whose columns represent the genes and rows represent patients
#### SRT data sets [Data SRT](SRT-data)
For each of the SRT samples, the spatially resolved transcriptomics (SRT) gene expression data matrix whose columns represent spots and rows represent genes, and spot locations denoting the spacial location of spots in the tissue slice

#### Results of the bulk data analysis [Result](result)
The following files proves the results for the analysis on bulk data of cancer patients.
* Genes are sorted by thier predictive values in the Eigen-Patients for the 5 phenotypes [Eigen-Patients](result/Eigen-Patients.xlsx).
* GO terms enriched for the sorted list of top predictor genes [GO terms for predictors](result/Enriched-GO-terms.xlsx). The GO terms enriched for top negative and positive predictors in Eigen-Patient.
#### Results of the SRT data analysis [Result](result)
The following files proves the results for the analysis on spacial transcriptomics (SRT) data
* Predicted phenotype quantity on spots [Predicted values](result/Prediction-CosineSimilarity.xlsx). Prediction based on cosine similarity between Eigen-Patient and SRT gene expression in each spot.
* Plots for showing the prediction results over spatial location [Prediction Plots](result)
