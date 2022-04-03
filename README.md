# BinaryClust
Tools for CyTOF analysis


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

