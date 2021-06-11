# Analyzing Gene Expression Data for  the cancer type Lung
Squamous Cell Carcinoma (LUSC)  

Data folder includes 2 files:

1. “lusc-rsem-fpkm-tcga-t_paired.txt”: GE data for tissues with cancer.
2. “lusc-rsem-fpkm-tcga_paired.txt” : GE data for tissues in a healthy case. 

- Each row in the two files represent a gene, and the columns represent the expression level of
  this gene in different samples.
- Data are paired meaning that a healthy sample and a diseased sample are taken from the same
  subject. So, both GE files have the same number of cases with the same order.
- Files are tab-separated.  

## Analysis steps:

### 1. Reading Files

```python
import pandas as pd
#----------------------------------------------reading the files----------------------------------------------#
healthy = pd.read_csv('data/lusc-rsem-fpkm-tcga_paired.txt', sep='\t')
cancer = pd.read_csv('data/lusc-rsem-fpkm-tcga-t_paired.txt', sep='\t')
pd.options.display.max_columns = None

# copying the data
h = healthy
c = cancer
```

### 2. Filteration 

```python
import numpy as np
#----------------------------------------------Files Filteration----------------------------------------------#

# drop the 2nd column that identify the genes
h.drop(["Entrez_Gene_Id"], axis = 1, inplace = True)
c.drop(["Entrez_Gene_Id"], axis = 1, inplace = True)

Clength = len(c.columns) # return number of columns after dropping the 2nd identification column
droplength = ( 0.5 * (Clength-1) +1) # return the number of 50% of total columns ,+1 because we still have the numeclature column

# replace the zeros (missing values) in file with nan this is required for next step
h = h.replace(0,np.nan)
c = c.replace(0,np.nan)

# in order to remove the rows that have more than 50% missing data dropna function is applied
# in which we put threshold that reflects the minimum number of not nan data required in each row
h = h.dropna(thresh = droplength) 
c = c.dropna(thresh = droplength)

#return back to zero values to proceed calculations
h = h.replace(np.nan,0)
c = c.replace(np.nan,0)

# Extracting the column of hugo_symbol from both files
G_h = h.iloc[0: , 0]
G_c = c.iloc[0: , 0]

# getting the intersected genes in both files according to hugo_symbol extracted from the previous step
common = pd.merge(G_h, G_c, how= 'inner')

# rewriting the 2 files with the common genes only so that each file has the no & name of genes  
h = pd.merge(h, common, how= 'inner', on= ['Hugo_Symbol'])
c = pd.merge(c, common, how= 'inner', on= ['Hugo_Symbol']) 

```



### 3. Linear Correlation (Pearson's Correlation) Calculation

```python
from scipy.stats import pearsonr
#--------------------------------pearson (linear) correlation calculation-------------------------------------#

data_corr = []  #a list saves index and correlation
corr_and_gene = []  #a list saves gene name and correlation
mycorr = []      #a list saves correlation value only
G_H = []        # a list saves the genes excluding the 1st column "gene_name" (healthy)
G_C = []        # a list saves the genes excluding the 1st column "gene_name" (cancer)

for x in range(0,len(h)): #len(h) = 17337
    G_h = h.iloc[x, 1:]
    G_c = c.iloc[x, 1:]
    G_H.append(G_h)
    G_C.append(G_c)
    name_gene = h.iloc[x, 0]
    r, _ = pearsonr(G_h, G_c)
    data_corr.append((x,r))
    corr_and_gene.append((name_gene,r))
    mycorr.append(r)
    
#get maximum corr and index   
max_value = max(mycorr) #Return the max value of the list.
max_index = mycorr. index(max_value) #Find the index of the max value.
print('max CC: ',max_value)
print('index of max CC: ',max_index)
max_gene_name = h.iloc[max_index, 0]
print("max correlation in gene",max_gene_name)
min_value = min(mycorr)
min_index = mycorr. index(min_value) #Find the index of the max value.
print('min CC: ',min_value)
print('index of min CC: ', min_index)
min_gene_name = h.iloc[min_index, 0]
print("min correlation in gene",min_gene_name)


#sorting the correlation values
mycorr.sort(reverse=True)
```



### 4. Plotting Correlation Results

```python
import matplotlib.pyplot as plt
#----------------------------------------------Plotting the results----------------------------------------------#

H_gh = h.iloc[max_index, 1:].astype(np.float)
H_gc = c.iloc[max_index, 1:].astype(np.float)
L_gh = h.iloc[min_index, 1:].astype(np.float)
L_gc = c.iloc[min_index, 1:].astype(np.float)

#first plot
plt.scatter(H_gh,H_gc, label='original data')
#curve fitting
m, b = np.polyfit(H_gh, H_gc, 1)
plt.plot(H_gh, m*H_gh + b, c = 'r', label='fitted line')
plt.title('The expression levels of genes which have max CC.')
plt.xlabel("Healthy")
plt.ylabel('Cancerous')
plt.legend()
plt.show()

#second plot
plt.scatter(L_gh,L_gc, label='original data')
#curve fitting
m, b = np.polyfit(L_gh, L_gc, 1)
plt.plot(L_gh, m*L_gh + b, c = 'r', label='fitted line')
plt.title('The expression levels of genes which have min CC.')
plt.xlabel("Healthy")
plt.ylabel('Cancerous')
plt.legend()
plt.show()
```



### 5. Hypothesis Testing  &  FDR Correction

```python
from scipy.stats import ttest_rel
from statsmodels.stats.multitest import multipletests
from scipy.stats import ttest_ind
#----------------------------------------------Hypothesis test----------------------------------------------# 

#Samples are Paired
P_val_paired = []  #list saves P values for paired samples 
for x in range(0,len(h)):
    Gh = G_H[x]
    Gc = G_C[x]
    p_val = ttest_rel(Gh, Gc).pvalue
    P_val_paired.append(p_val)

#Multiple tests correction (FDR Correction) for paired samples
corrected_P_val_paired = multipletests(P_val_paired, alpha=0.05, method='fdr_bh')[1]

#Samples are independent 
P_val_independent = []  #list saves P values for paired samples 
for x in range(0,len(h)):
    Gh = G_H[x]
    Gc = G_C[x]
    p_val = ttest_ind(Gh,Gc).pvalue
    P_val_independent.append(p_val)

#Multiple tests correction (FDR Correction) for independent samples
corrected_P_val_independent = multipletests(P_val_independent, alpha=0.05, method='fdr_bh')[1]

# significance_genes paired
gene_names = h.iloc[0: , 0]
significance_genes_paired = pd.DataFrame({'Gene_name':gene_names, 'p-values':P_val_paired, 'p-values_fdr':corrected_P_val_paired})

significance_genes_paired['significance:p_values'] = significance_genes_paired['p-values'].apply(lambda x: x < 0.05)
significance_genes_paired['significance:p_values_fdr'] = significance_genes_paired['p-values_fdr'].apply(lambda x: x < 0.05)

# Get significant genes before fdr correction
diffrentially_genes_paired_b = significance_genes_paired[significance_genes_paired['significance:p_values']== True]

# To get gene names before fdr
diff_paired_b=diffrentially_genes_paired_b['Gene_name'].to_list()

#convert this list to dataframe and save it in a file
# df_diff_paired_b = pd.DataFrame(diff_paired_b)
# df_diff_paired_b.to_csv(r'C:/Users/user/Desktop\DEGS_paired_before_FDR.csv')

#Get significant genes after fdr correction
diffrentially_genes_paired = significance_genes_paired[significance_genes_paired['significance:p_values_fdr']== True]

# To get gene names after fdr
diff_paired=diffrentially_genes_paired['Gene_name'].to_list()

#convert this list to dataframe and save it in a file
# df_diff_paired = pd.DataFrame(diff_paired)
# df_diff_paired.to_csv(r'C:/Users/user/Desktop\DEGS_paired_after_FDR.csv')


# significance_genes independent
significance_genes_independent = pd.DataFrame({'Gene_name':gene_names, 'p-values':P_val_independent, 'p-values_fdr':corrected_P_val_independent})
significance_genes_independent['significance:p_values'] = significance_genes_independent['p-values'].apply(lambda x: x < 0.05)
significance_genes_independent['significance:p_values_fdr'] = significance_genes_independent['p-values_fdr'].apply(lambda x: x < 0.05)

# Get significant genes before fdr correction
diffrentially_genes_independent_b = significance_genes_independent[significance_genes_independent['significance:p_values']== True]

# To get gene names before fdr
diff_independent_b=diffrentially_genes_independent_b['Gene_name'].to_list()

#convert this list to dataframe and save it in a file
# df_diff_independent_b = pd.DataFrame(diff_independent_b)
# df_diff_independent_b.to_csv(r'C:/Users/user/Desktop\DEGS_independent_before_FDR.csv')

#Get significant genes after fdr correction
diffrentially_genes_independent = significance_genes_independent[significance_genes_independent['significance:p_values_fdr']== True]

# To get gene names after fdr
diff_independent=diffrentially_genes_independent['Gene_name'].to_list()

#convert this list to dataframe and save it in a file
# df_diff_independent = pd.DataFrame(diff_independent)
# df_diff_independent.to_csv(r'C:/Users/user/Desktop\DEGS_independent_after_FDR.csv')


#Get The distinct values in diff_paired (Not in diff_independent ) after FDR
diff_paired_distinct = [] #a list saves distinct values in DEGs paired sets
DEGs_common = []  #a list saves common values between  DEGs paired sets and independent paired sets
for x in diff_paired:
  if x not in diff_independent:
    diff_paired_distinct.append(x)
  elif x  in diff_independent:
    DEGs_common.append(x)

#Get The distinct values in diff_independent (Not in diff_paired) after FDR
diff_independent_distinct = [] #a list saves distinct values in DEGs independent sets
for x in diff_independent:
  if x not in diff_paired:
    diff_independent_distinct.append(x)

#Making sure that no common elements between diff_paired_distinct and diff_independent
compare = [] #a list saves an element count in diff_independent (if exists)
for i in diff_paired_distinct:
   compare.append(diff_independent.count(i))

#Making sure that no common elements between  diff_independent_distinct and diff_paired
Compare = [] #a list saves an element count in diff_paired (if exists)
for i in diff_independent_distinct:
   Compare.append(diff_paired.count(i))


if len(DEGs_common)+len(diff_paired_distinct)==len(diff_paired) and len(DEGs_common)+len(diff_independent_distinct)==len(diff_independent) :
    print('The Comparison has been done successfuly' )
    
```

#### Saving Outputs in to SpreadSheets

```python
#save diff_independent_distinct ,diff_paired_distinct ,DEGs_common in spreadsheet

df_diff_paired_distinct = pd.DataFrame(diff_paired_distinct)
df_diff_paired_distinct.to_csv(r'C:/Users/user/Desktop\DEGS_paired_distinct_after_FDR.csv')

df_diff_independent_distinct = pd.DataFrame(diff_independent_distinct)
df_diff_independent_distinct.to_csv(r'C:/Users/user/Desktop\DEGS_independent_distinct_after_FDR.csv')

df_DEGs_common = pd.DataFrame(DEGs_common)
df_DEGs_common.to_csv(r'C:/Users/user/Desktop\DEGS_common_after_FDR.csv')
```

