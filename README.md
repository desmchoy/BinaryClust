# BinaryClust
Tools for CyTOF analysis



## Introduction
CyTOF data is log normal with zero inflation for most markers. For example:

![CD3](/images/CD3.png)

BinaryClust takes advantage of this feature to subset cells for clustering by performing binary classifications on the markers. Using _K_-means clustering (_K_ = 2), it is possible to automatically separate cells into positive and negative populations. Here, _K_-means breaks down the cells into CD3+ (red) and CD3- (blue) populations:

![CD3 K-means](/images/CD3_shaded.png)

This process is done in BinaryClust on all the markers separately to generate a classification matrix. A user-defined selection criteria will then classify these cells into different populations, which will then be individually clustered:

![BinaryClust Schematic](/images/schematic.png)

BinaryClust also comes with a set of functions to perform standard analysis tasks such as data exploration and data visualisation by _t_-SNE or UMAP.





## Prerequisites and Installation
To intall BinaryClust, you will first need to install the R packages devtools and BiocManager:

```
install.packages('devtools')
install.packages('BiocManager')
```

Three Bioconductor packages are needed:

```
library(BiocManager)
install('ComplexHeatmap')
install('flowCore')
install('FlowSOM')
```

Then you can install BinaryClust and its dependencies:

```
library(devtools)
install_github("desmchoy/BinaryClust")
```

Don't forget to load the package after installation:

```
library(BinaryClust)
```

The original BinaryClust pipeline was designed for large-scale parallel runs on a high performance computing (HPC) facility. However, for reasonably sized samples, it can be used on local installations or via an IDE such as RStudio. though it is not recommended to run _t_-SNE or UMAP for large samples on a local installation due to computational demands.





## Quick Start
BinaryClust is designed to be efficient and straight to the point. A quick analysis can be run within seconds using a single command given an FCS file and a user-defined cell classification file in CSV format (see the following section for details):

```
BC.results <- run_BinaryClust('HD.fcs', 'cell_types.csv')
```

`run_BinaryClust` is a shortcut function that performs five individual steps of analysis (each elaborated in the following section) using default settings in one go and returns their results in one single object. The result is a list of four data frames that provides the inputs for plotting functions or further analysis. The first data frame gives you arcsinh transformed (cofactor 5) data:

```
> head(BC.results[[1]])
      CD45 CD196_CCR6      CD19 CD127_IL-7Ra       CD38       CD33        IgD
1 4.763022 0.00000000 0.0000000   0.00000000 0.02726449 0.00000000 0.00000000
2 3.912244 2.05420158 0.0000000   2.45454068 0.01258484 0.00000000 0.03181502
3 5.385946 0.00000000 0.0000000   1.71571358 0.00000000 0.21548775 0.00000000
4 4.951837 0.44154719 0.0000000   0.06832929 0.26353484 0.00000000 0.00000000
5 5.501669 2.19138813 0.8805520   0.31980510 3.11632848 2.74392043 0.12969393
6 5.213904 0.09108109 0.1907823   2.25250901 1.42557094 0.06746715 0.20465348
      CD11c      CD16 CD194_CCR4     CD34 CD123_IL-3R     TCRgd CD185_CXCR5
1 0.0000000 0.0000000  0.0000000 0.000000     0.00000 0.0000000           0
2 0.3241666 0.2583896  2.5192846 0.000000     0.00000 0.0000000           0
3 0.0000000 0.0000000  0.0000000 0.000000     0.11153 0.0000000           0
4 0.1237903 0.0000000  0.0000000 0.000000     0.00000 0.7061377           0
5 5.0112698 0.4296513  0.8198227 1.515479     0.00000 0.1518769           0
6 0.0000000 0.0000000  0.4235676 0.000000     0.00000 0.3448717           0
...
```

The second data frame gives you the median expressions of each marker and the abundance of each cell type definied in the cell classification file:

```
> head(BC.results[[2]])
        Cell.Type     CD45 CD196_CCR6       CD19 CD127_IL-7Ra      CD38
1         B Cells 5.175686  2.7822063 4.55294083   0.24782457 2.6548137
2 Dendritic Cells 4.973808  0.1735834 0.05988902   0.00000000 2.8594343
3       Monocytes 5.068816  0.2541294 0.07413499   0.06794235 3.6651274
4        NK Cells 4.780354  0.0000000 0.05407581   0.00000000 3.9153949
5    T Cells, CD4 5.146704  0.2767183 0.00000000   2.36491443 0.5727581
6    T Cells, CD8 5.126842  0.0000000 0.00000000   0.37595694 0.2547606
       CD33       IgD      CD11c       CD16   CD194_CCR4        CD34
1 0.1094124 4.3819244 0.08709911 0.02833477 0.0000000000 0.000000000
2 1.6333382 0.1881026 4.83674662 1.86184415 0.3341434080 0.212236144
3 2.4676174 0.2542666 4.75456686 1.06149001 0.4275743442 0.213740130
4 0.1222321 0.1060854 3.17299179 3.55946059 0.0002289466 0.007264909
5 0.0000000 0.0000000 0.00000000 0.00000000 0.5755964856 0.059087905
6 0.0000000 0.0000000 0.02973928 0.00000000 0.0000000000 0.000000000
...
CD4      CD14  CD56_NCAM       CD11b Frequency Percentage
1 0.8065952 0.0000000 0.00000000 0.000000000      4105   6.400262
2 3.2867254 0.4042923 0.07917726 0.054064754      7765  12.106707
3 3.5688402 2.8055202 0.19211170 0.058509953     18883  29.441205
4 0.1929261 0.0000000 4.05164056 0.007594521      5561   8.670367
5 5.4276076 0.3077517 0.25694324 0.000000000     11667  18.190464
6 0.1380546 0.0000000 0.16847492 0.000000000     14128  22.027503

```

The third data frame lists out exhaustively the clustering results of each cell type:

```
> head(BC.results[[3]])
     Cell.Type Cluster
1 T Cells, CD8       8
2 T Cells, CD4      32
3 T Cells, CD8       4
4 T Cells, CD8       9
5    Monocytes      30
6 T Cells, CD4      19
...
```

The fourth data frame summarises the clustering results. It gives the median expressions of each marker of each cluster of each cell type defined in the cell classification file and their abundances. It also leaves an empty `Cell.Subtype` column for filling for further analysis.

```
> head(BC.results[[4]])
  Cell.Type Cluster     CD45 CD196_CCR6     CD19 CD127_IL-7Ra      CD38
1   B Cells       1 5.169841   2.830101 4.535509   0.33521943 3.3112159
2   B Cells      10 5.088842   2.498848 4.502411   0.17985776 2.3036351
3   B Cells      11 5.400011   2.562452 4.672512   0.19210797 0.9343564
4   B Cells      12 5.389607   5.242463 4.983812   1.87714720 3.1863538
5   B Cells      13 5.360868   2.147824 4.561302   0.11040360 0.6605718
6   B Cells      14 4.105287   2.512916 4.147211   0.04629935 0.9779664
         CD33        IgD      CD11c      CD16 CD194_CCR4        CD34
1 0.289654941 5.41413151 0.09756387 0.1781481          0 0.003571104
2 0.000000000 0.00000000 0.04685361 0.0000000          0 0.000000000
3 0.006471317 3.81854993 0.00000000 0.0000000          0 0.000000000
4 0.238622565 5.00851295 0.07870780 0.1452360          0 0.020131121
5 0.000000000 0.07353001 0.03799469 0.0000000          0 0.000000000
6 0.000000000 2.67294214 0.00000000 0.0000000          0 0.000000000
...
CD14   CD56_NCAM       CD11b Frequency Percentage Cell.Subtype
1 0.000000000 0.006773511 0.000000000       164 0.25569865
2 0.000000000 0.000000000 0.000000000       135 0.21048364
3 0.000000000 0.000000000 0.000000000       134 0.20892451
4 0.000000000 0.000000000 0.000000000        37 0.05768811
5 0.004020624 0.000000000 0.005794109       114 0.17774174
6 0.000000000 0.000000000 0.006619664        79 0.12317191
```

These data frames can be used to generate plots using the relevant functions of `BinaryClust`. For example, the binary classification summary (second data frame) can be plotted out:

```
plot_binary_abundances(BC.results[[2]])
plot_binary_median(BC.results[[2]])
```

Each of these analysis steps will now be explained in the following section.



## Data Exploration



## Required Input Files
To run standard BinaryClust analysis for a single CyTOF file, two files are needed:

1. A CyTOF FCS file (formats prior to FCS 3.1 not tested)
2. A cell classification file in CSV format

##### The Cell Classification File
The cell stratification file should be in CSV format. An example is shown below. 
```
Cell Type,CD14,CD16,CD161,CD19,CD20,CD3,CD4,CD56_NCAM,CD8,TCRgd
NK Cells,-,A,A,-,A,-,A,+,A,A
Dendritic Cells,-,A,A,-,A,-,A,-,A,A
Monocytes,+,A,A,-,-,-,A,A,A,A
B Cells,-,-,-,+,A,-,A,A,A,A
"T Cells, Gamma Delta",-,A,A,-,A,+,A,A,A,+
"T Cells, CD4",-,A,A,-,A,+,+,A,-,-
"T Cells, CD8",-,A,A,-,A,+,-,A,+,-
```





