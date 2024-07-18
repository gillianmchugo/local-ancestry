# bite bootstrap
## generate treemix bootstrap script with bite package
### install required packages
if (!require(beanplot)) install.packages("beanplot")
if (!require(car)) install.packages("car")
if (!require(devtools)) install.packages("devtools")
if (!require(kableExtra)) install.packages("kableExtra")
if (!require(knitr)) install.packages("knitr")
if (!require(RColorBrewer)) install.packages("RColorBrewer")
if (!require(reshape2)) install.packages("reshape2")
if (!require(rmarkdown)) install.packages("rmarkdown")
devtools::install_url("https://cran.r-project.org/src/contrib/Archive/GenABEL.data/GenABEL.data_1.0.0.tar.gz")
devtools::install_url("https://cran.r-project.org/src/contrib/Archive/GenABEL/GenABEL_1.8-0.tar.gz")
install.packages("https://cran.r-project.org/src/contrib/Archive/RCircos/RCircos_1.1.3.tar.gz",
                 repos = NULL,
                 method = "wget")
install.packages("/path/to/bite/BITE-master/BITE_1.2.0008.tar.gz",
                 repos = NULL,
                 type = "source")

### load bite package
library(BITE)

### get current working directory
home <- getwd()

### set target working directory
setwd("/home/workspace/gmchugo/local-ancestry/06_phylogenetic-analysis/01_input")

### generate treemix bootstrap script
treemix.scripts()

### reset working directory to current directory
setwd(home)