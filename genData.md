# Generate test data

This file actually generates the data used in the examples, and for testing. It should be re-run every time packages are updated so that comparison of results are valid. This will be incredibly important with new versions of Bioconductor.

## Set up knitr and R options


```r
knit_hooks$set(inline = identity)  # sets output for inline resutls nicely
knit_hooks$set(error = function(x, options) stop(x))  # kills knitr if there
# is an error, therefore we don't waste time generating error messages
options(save.defaults = list(compress = "xz"), stringsAsFactors = FALSE)
```




## Load libraries


```r
library("affy")
```

```
## Loading required package: BiocGenerics
```

```
## Attaching package: 'BiocGenerics'
```

```
## The following object(s) are masked from 'package:stats':
## 
## xtabs
```

```
## The following object(s) are masked from 'package:base':
## 
## anyDuplicated, cbind, colnames, duplicated, eval, Filter, Find, get,
## intersect, lapply, Map, mapply, mget, order, paste, pmax, pmax.int, pmin,
## pmin.int, Position, rbind, Reduce, rep.int, rownames, sapply, setdiff,
## table, tapply, union, unique
```

```
## Loading required package: Biobase
```

```
## Welcome to Bioconductor
## 
## Vignettes contain introductory material; view with 'browseVignettes()'. To
## cite Bioconductor, see 'citation("Biobase")', and for packages
## 'citation("pkgname")'.
```

```r
library("hgu95av2.db")
```

```
## Loading required package: AnnotationDbi
```

```
## Loading required package: org.Hs.eg.db
```

```
## Loading required package: DBI
```

```
## NA
```

```
## NA
```

```r
library("hgu95av2cdf")
```

```
## NA
```

```r
library("genefilter")
library("estrogen")
library("limma")
library("categoryCompare")
```

```
## Loading required package: Category
```

```
## NA
```

```r
library("GO.db")
```






```r
# read in the raw estrogen data
datadir <- system.file("extdata", package = "estrogen")
pd <- read.AnnotatedDataFrame(file.path(datadir, "estrogen.txt"), 
    header = TRUE, sep = "", row.names = 1)
pData(pd)
```

```
##              estrogen time.h
## low10-1.cel    absent     10
## low10-2.cel    absent     10
## high10-1.cel  present     10
## high10-2.cel  present     10
## low48-1.cel    absent     48
## low48-2.cel    absent     48
## high48-1.cel  present     48
## high48-2.cel  present     48
```

```r

a <- ReadAffy(filenames = file.path(datadir, rownames(pData(pd))), 
    phenoData = pd, verbose = TRUE)
```

```
## 1 reading J:/R150_libraries/estrogen/extdata/low10-1.cel ...instantiating an AffyBatch (intensity a 409600x8 matrix)...done.
## Reading in : J:/R150_libraries/estrogen/extdata/low10-1.cel
## Reading in : J:/R150_libraries/estrogen/extdata/low10-2.cel
## Reading in : J:/R150_libraries/estrogen/extdata/high10-1.cel
## Reading in : J:/R150_libraries/estrogen/extdata/high10-2.cel
## Reading in : J:/R150_libraries/estrogen/extdata/low48-1.cel
## Reading in : J:/R150_libraries/estrogen/extdata/low48-2.cel
## Reading in : J:/R150_libraries/estrogen/extdata/high48-1.cel
## Reading in : J:/R150_libraries/estrogen/extdata/high48-2.cel
```

```r
eData <- rma(a)
```

```
## Background correcting
## Normalizing
## Calculating Expression
```




## Set up time point 10


```r
pCut <- 0.05
e10 <- eData[, eData$time.h == 10]
e10 <- nsFilter(e10, remove.dupEntrez = TRUE, var.filter = FALSE, 
    feature.exclude = "^AFFX")$eset

e10$estrogen <- factor(e10$estrogen)
d10 <- model.matrix(~0 + e10$estrogen)
colnames(d10) <- unique(e10$estrogen)
fit10 <- lmFit(e10, d10)
c10 <- makeContrasts(present - absent, levels = d10)
fit10_2 <- contrasts.fit(fit10, c10)
eB10 <- eBayes(fit10_2)
table10 <- topTable(eB10, number = nrow(e10), p.value = pCut, adjust.method = "BH")
table10$Entrez <- unlist(mget(table10$ID, hgu95av2ENTREZID, ifnotfound = NA))
```




## Set up time point 48


```r
e48 <- eData[, eData$time.h == 48]
e48 <- nsFilter(e48, remove.dupEntrez = TRUE, var.filter = FALSE, 
    feature.exclude = "^AFFX")$eset

e48$estrogen <- factor(e48$estrogen)
d48 <- model.matrix(~0 + e48$estrogen)
colnames(d48) <- unique(e48$estrogen)
fit48 <- lmFit(e48, d48)
c48 <- makeContrasts(present - absent, levels = d48)
fit48_2 <- contrasts.fit(fit48, c48)
eB48 <- eBayes(fit48_2)
table48 <- topTable(eB48, number = nrow(e48), p.value = pCut, adjust.method = "BH")
table48$Entrez <- unlist(mget(table48$ID, hgu95av2ENTREZID, ifnotfound = NA))
```




## Gene Universe


```r
gUniverse <- unique(unlist(mget(rownames(as.matrix(exprs(e48))), 
    hgu95av2ENTREZID, ifnotfound = NA)))
gUniverse <- gUniverse[!(is.na(gUniverse))]
```




# Basic calculations

Perform some of the basic calculations that are needed so we can make a comparison between the results of test operations and "truth", *i.e.* so that we know what is working. Note, for **ccData**, we only keep the top 100 genes in **T10** and **T48**, this keeps the filesize down and the required time for calculations.



```r
table48 <- table48[1:100, ]
table10 <- table10[1:100, ]

g10 <- unique(table10$Entrez)
g10 <- g10[!(is.na(g10))]
g48 <- unique(table48$Entrez)
g48 <- g48[!(is.na(g48))]

hyperGParamT <- new("GOHyperGParamsCC", geneIds = g10, universeGeneIds = gUniverse, 
    annotation = "org.Hs.eg.db", ontology = "CC", conditional = FALSE, testDirection = "over", 
    fdr = 0, pvalueCutoff = 0.05)

list10 <- list(genes = g10, universe = gUniverse, annotation = "org.Hs.eg.db")
list48 <- list(genes = g48, universe = gUniverse, annotation = "org.Hs.eg.db")

geneLists <- list(T10 = list10, T48 = list48)
geneLists <- new("ccGeneList", geneLists, ccType = "CC", fdr = 0)
```

```
## Warning: NAs introduced by coercion
```

```r

enrichLists <- ccEnrich(geneLists)
```

```
## Performing Enrichment Calculations ....
## T10 : CC 
## T48 : CC 
## Done!!
## 
```

```r
ccOpts <- new("ccOptions", listNames = names(geneLists))
ccResults <- ccCompare(enrichLists, ccOpts)
graphType(enrichLists$CC) <- "hierarchical"
ccResultsHier <- ccCompare(enrichLists$CC, ccOpts)
```




# Save data


```r
.sessionInfo <- sessionInfo()
.timeDate <- Sys.time()

saveList <- c("table48", "table10", "gUniverse", "hyperGParamT", 
    "geneLists", "enrichLists", "ccOpts", "ccResults", "ccResultsHier", ".sessionInfo", 
    ".timeDate")

save(list = saveList, file = "data/ccData.RData")
cat("ccData:", saveList, file = "data/datalist", sep = " ")
.sessionInfo
```

```
## R version 2.15.0 (2012-03-30)
## Platform: x86_64-pc-mingw32/x64 (64-bit)
## 
## locale:
## [1] LC_COLLATE=English_United States.1252 
## [2] LC_CTYPE=English_United States.1252   
## [3] LC_MONETARY=English_United States.1252
## [4] LC_NUMERIC=C                          
## [5] LC_TIME=English_United States.1252    
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] GO.db_2.7.1           categoryCompare_1.1.0 Category_2.22.0      
##  [4] limma_3.12.1          estrogen_1.8.8        genefilter_1.38.0    
##  [7] hgu95av2cdf_2.10.0    hgu95av2.db_2.7.1     org.Hs.eg.db_2.7.1   
## [10] RSQLite_0.11.1        DBI_0.2-5             AnnotationDbi_1.18.1 
## [13] affy_1.34.0           Biobase_2.16.0        BiocGenerics_0.2.0   
## [16] knitr_0.6            
## 
## loaded via a namespace (and not attached):
##  [1] affyio_1.24.0         annotate_1.34.0       BiocInstaller_1.4.6  
##  [4] codetools_0.2-8       colorspace_1.1-1      digest_0.5.2         
##  [7] evaluate_0.4.2        formatR_0.4           GOstats_2.22.0       
## [10] graph_1.34.0          GSEABase_1.18.0       highlight_0.3.1      
## [13] hwriter_1.3           IRanges_1.14.3        parser_0.0-14        
## [16] plyr_1.7.1            preprocessCore_1.18.0 RBGL_1.32.0          
## [19] Rcpp_0.9.10           RCurl_1.91-1.1        RCytoscape_1.6.3     
## [22] splines_2.15.0        stats4_2.15.0         stringr_0.6          
## [25] survival_2.36-14      tools_2.15.0          XML_3.9-4.1          
## [28] XMLRPC_0.2-4          xtable_1.7-0          zlibbioc_1.2.0       
```

```r
.timeDate
```

```
## [1] "2012-06-14 13:39:19 EDT"
```



