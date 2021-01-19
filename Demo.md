---
title: <font size = "6">Demo of microbial fcm analysis pipeline</font>
output:
  html_document:
    code_folding: show
    highlight: haddock
    keep_md: yes
    theme: flatly
    toc: yes
    number_sections: true
    toc_float:
      collapsed: no
      smooth_scroll: yes
      toc_depth: 3
editor_options:
  chunk_output_type: console
---

<style type="text/css">
h1.title {
  font-size: 38px;
  color: black;
  text-align: center;
}

</style>



# Libraries


```r
library("Phenoflow")
library("plyr")
library("dplyr")
library("ggplot2")
library("flowAI")
library("scales")
library("plyr")
library("cowplot")
library("ggcyto")
library("RColorBrewer")
library("grid")
library("tidyr")
library("flowFP")
library("FlowSOM")
library("readxl")
# for flowEMMI
library("Rcpp")
library("RcppEigen")
library("mixtools")
library("gtools")
library("flowCore")
library("flowViz")
library("randomcoloR")
sourceCpp("ext/flowEMMi.cpp")
```

```
## 
## > flowEMMi <- function(frame, ch1 = "FS.Log", ch2 = "FL.4.Log", 
## +     use_log = TRUE, diff.ll = 1, sample_size = 10, start_cluster = 8, 
## +     end_cl .... [TRUNCATED]
```

```r
# Set seed for reproducible analysis
set.seed(777)

# For avoiding verbose console output
run_quiet <- function (x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

# Import metadata
metadata <- readxl::read_excel("./metadata.xlsx")
metadata$Sample_names <- gsub(".5", "_5", metadata$Sample_names, fixed = TRUE)
metadata$Sample_names <- gsub(".fcs", "_QC.fcs", metadata$Sample_names, fixed = TRUE)
metadata$Treatment <- factor(metadata$Treatment,
                             levels = c("Adaptation of AMC",
                                        "AMC after temperature disturbance",
                                        "AMC after pH disturbance",
                                        "Adaptation of CMC",
                                        "CMC after temperature disturbance",
                                        "CMC after pH disturbance"))
```

# Preprocessing

## Import data 

Data set from [Liu et al (2017).](https://msphere.asm.org/content/3/1/e00564-17)


```r
# Import data
flowData <- read.flowSet(path = "data/zishu_2018/")
param <- c("FS Log", "SS Log","FL 4 Log")
sampleNames(flowData) <- gsub(".5", "_5", sampleNames(flowData), fixed = TRUE)
```

## Transform data 


```r
# Transform data
flowData_transformed_log <- transform(flowData,`FS Log`=log10(`FS Log`),
                                   `SS Log`=log10(`SS Log`),
                                   `FL 4 Log`=log10(`FL 4 Log`))
flowData_transformed_asinh <- transform(flowData,`FS Log`=asinh(`FS Log`),
                                   `SS Log`=asinh(`SS Log`),
                                   `FL 4 Log`=asinh(`FL 4 Log`))
```

## Denoise data - step 1 {.tabset}

### Untransformed data 


```r
r_sam <- sample(1:length(flowData), 6)

#  scatter plot
p_scatter1 <- ggcyto::ggcyto(flowData[r_sam], aes(x = `FS Log`, y = `FL 4 Log`)) + 
  geom_hex(bins = 300) +
  theme_bw()+
  labs(x = "Forward scatter (a.u.)", y = "DAPI (a.u.)")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "white"))
  # coord_cartesian(xlim = c(0,6), ylim = c(0,6))

print(p_scatter1)
```

<img src="./Figures/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />

### Log-transformed data 


```r
#  scatter plot
p_scatter2 <- ggcyto::ggcyto(flowData_transformed_log[r_sam], aes(x = `FS Log`, y = `FL 4 Log`)) + 
  geom_hex(bins = 300) +
  theme_bw()+
  labs(x = "Forward scatter (a.u.)", y = "DAPI (a.u.)")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "white"))
  # coord_cartesian(xlim = c(0,6), ylim = c(0,6))

print(p_scatter2)
```

<img src="./Figures/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

### asinh-transformed data 


```r
#  scatter plot
p_scatter3 <- ggcyto::ggcyto(flowData_transformed_asinh[r_sam], aes(x = `FS Log`, y = `FL 4 Log`)) + 
  geom_hex(bins = 300) +
  theme_bw()+
  labs(x = "Forward scatter (a.u.)", y = "DAPI (a.u.)")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "white"))
  # coord_cartesian(xlim = c(0,6), ylim = c(0,6))

print(p_scatter3)
```

<img src="./Figures/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

## Denoise data - step 2 {.tabset}

### Without gate


```r
#  scatter plot
p_scatter4 <- ggcyto::ggcyto(flowData_transformed_asinh[r_sam], aes(x = `FS Log`, y = `FL 4 Log`)) + 
  geom_hex(bins = 300) +
  theme_bw()+
  labs(x = "Forward scatter (a.u.)", y = "DAPI (a.u.)")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "white"))
  # coord_cartesian(xlim = c(0,6), ylim = c(0,6))

print(p_scatter4)
```

<img src="./Figures/unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

### With gate


```r
# Create a PolygonGate for denoising the dataset
polycut <- matrix(c(1.5, 1.5, 4, 4, 9, 9, 6.75, 6.75, 6.25, 6.25,
  2.25, 6, 6, 8.75, 8.75, 2.25, 2.25, 2.75, 2.75, 2.25),
     ncol = 2,
     nrow = 10)
colnames(polycut) <- c("FS Log", "FL 4 Log")
polyGate <- polygonGate(.gate = polycut, filterId = "Cell population")

#  scatter plot
p_scatter5 <- ggcyto::ggcyto(flowData_transformed_asinh[r_sam], aes(x = `FS Log`, y = `FL 4 Log`)) + 
  geom_hex(bins = 300) +
  theme_bw()+
  labs(x = "Forward scatter (a.u.)", y = "DAPI (a.u.)")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "white"))+
  # coord_cartesian(xlim = c(0,6), ylim = c(0,6))+
  geom_gate(polyGate, col = "#CB0001", fill = "#ffa401", alpha = 0.8, size = 1)
  
print(p_scatter5)
```

<img src="./Figures/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

```r
# Retain only events in gate
flowData_transformed_asinh <- Subset(flowData_transformed_asinh, polyGate)
```

## Denoise data - step 3

Doublet/clump removal not possible for this dataset due to lack of combined
-A/-H parameters.

## Denoise data - step 4 {.tabset}

### Before flowAI


```r
#  scatter plot
p_scatter6 <- ggcyto::ggcyto(flowData_transformed_asinh[r_sam], aes(x = `TIME`, y = `FL 4 Log`)) + 
  geom_hex(bins = 300) +
  theme_bw()+
  labs(x = "Time (ms)", y = "DAPI (a.u.)")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "white"))

print(p_scatter6)
```

<img src="./Figures/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />



### After flowAI 


```r
#  scatter plot
p_scatter7 <- ggcyto::ggcyto(flowData_transformed_asinh_dn[r_sam], aes(x = `TIME`, y = `FL 4 Log`)) + 
  geom_hex(bins = 300) +
  theme_bw()+
  labs(x = "Time (ms)", y = "DAPI (a.u.)")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "white"))

print(p_scatter7)
```

<img src="./Figures/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

### QC stats


```r
QC_stats <- read.delim("./data/QC_flowAI/QCmini.txt")

p_qc1 <- QC_stats %>%
  arrange(-X..anomalies) %>% 
  mutate(Name.file = factor(Name.file, levels = unique(Name.file))) %>% 
  ggplot(., aes(y = X..anomalies, x = Name.file))+
  geom_point(shape = 21, size = 3, fill = "#ffa401", alpha = 0.5)+
  theme_bw()+
  theme(axis.text = element_text(size = 12),
    axis.text.x = element_blank(),
        axis.title = element_text(size = 12))+
  labs(y = "% anomalies",  x = "Ranked samples")+
  geom_hline(yintercept = 5, linetype = 2)

p_qc2 <- QC_stats %>% 
  ggplot(., aes(x = X..anomalies))+
  geom_histogram(fill = "#ffa401", alpha = 0.5, col = "black")+
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
  labs(x = "% anomalies",  y = "Number of samples")+
  geom_vline(xintercept = median(QC_stats$X..anomalies), linetype = 2)+
  scale_x_continuous(breaks = seq(from = 0, to = 100, by = 10), limits = c(-5,100))

cowplot::plot_grid(p_qc1, p_qc2, align = "hv", ncol = 1)
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

```
## Warning: Removed 2 rows containing missing values (geom_bar).
```

<img src="./Figures/unnamed-chunk-12-1.png" style="display: block; margin: auto;" />

# Cell density measurement

Cell concentration estimates cannot be easily estimated for this dataset due to 
a lack of volumetric measurements. We refer the reader to the original manuscript.



# Cell population identification {.tabset}

## FlowSOM - fig1

```r
fSOM <- FlowSOM(
  input = flowData_transformed_asinh_dn[1:3],
  # Input options:
  compensate = FALSE,
  transform = FALSE,
  scale = FALSE,
  # SOM options:
  colsToUse = param,
  xdim = 7,
  ydim = 7,
  nClus = 10
  )

PlotStars(fSOM[[1]], backgroundValues = as.factor(fSOM[[2]]))
```

<img src="./Figures/unnamed-chunk-14-1.png" style="display: block; margin: auto;" />

## FlowSOM - fig2

```r
PlotMarker(fSOM[[1]], param[3])
```

<img src="./Figures/unnamed-chunk-15-1.png" style="display: block; margin: auto;" />

# Cytometric fingerprinting {.tabset}

## PhenoGMM


```r
fp_gmm <-
  PhenoGMM(
    flowData_transformed_asinh_dn,
    param = param,
    downsample = 1e2,
    auto_nG = TRUE,
    nG = 32,
    diagnostic_plot = TRUE
  )
```

```
## Your samples range between 235367 and 249862 cells
## Your samples were randomly subsampled to 100 cells
```

```r
# Calculate diversity indices
div_gmm <- Diversity_gmm(fp_gmm, R = 10)
```

```
## Thu Sep 03 12:44:17 2020 	Done with all 140 samples
```

```r
# Merge with metadata
div_gmm <- left_join(div_gmm, metadata, by = "Sample_names")
```


```r
# Plot results
div_gmm %>% 
  ggplot(., aes(x = Timepoint,  y = D2, fill = Treatment))+
  geom_point(shape = 21, size = 4)+
  geom_line()+
  facet_wrap(.~Treatment, scales = "free_x", ncol = 3)+
  scale_fill_brewer("", palette = "Accent")+
  theme_bw()+
  theme(strip.text = element_text(size = 8))+
  labs(y = "Hill diversity index (D2)")
```

<img src="./Figures/unnamed-chunk-17-1.png" style="display: block; margin: auto;" />

## PhenoFlow


```r
# Summary of max intensity values across data set
summary <-
  fsApply(
    x = flowData_transformed_asinh_dn,
    FUN = function(x)
      apply(x, 2, max),
    use.exprs = TRUE
  )
maxval <- max(summary[, "FL 4 Log"])

# Normalize intensities to [0,1] range for density estimation
mytrans <- function(x)
  x / maxval
flowData_transformed_asinh_dn_trans <-
  transform(
    flowData_transformed_asinh_dn,
    `FS Log` = mytrans(`FS Log`),
    `SS Log` = mytrans(`SS Log`),
    `FL 4 Log` = mytrans(`FL 4 Log`)
  )

# Estimate kernel densities across grids
fp_grid <-
  flowBasis(
    FCS_resample(flowData_transformed_asinh_dn_trans, 1e3),
    param = param,
    nbin = 128
  )
```

```
## Your samples range between 235367 and 249862 cells
## Your samples were randomly subsampled to 1000 cells
```

```r
# Perform ordination analysis
beta_fp <- beta_div_fcm(x = fp_grid, ord.type = "PCoA")
beta_fp <- data.frame(Sample_names = rownames(beta_fp$points), beta_fp$points)
beta_fp <- left_join(beta_fp, metadata, by = "Sample_names")
```


```r
# Plot results
beta_fp %>% 
  ggplot(., aes(x = X1,  y = X2, fill = Treatment))+
  geom_point(shape = 21, size = 4)+
  # facet_wrap(.~Treatment, ncol = 3)+
  scale_fill_brewer("", palette = "Accent")+
  theme_bw()+
  labs(x = "PCoA 1", y = "PCoA 2")+
  theme(strip.text = element_text(size = 8))
```

<img src="./Figures/unnamed-chunk-19-1.png" style="display: block; margin: auto;" />

## FlowFP

FlowFP does not seem to run in the current R and RStudio version. Example code
below in any case.


```r
# fp_FP <-
#   flowFPModel(
#     fcs = FCS_resample(flowData_transformed_asinh_dn, 1e3),
#     param = param,
#     nRecursions = 8
#   )
```

## FlowEMMi


```r
# fp_emmi <-
#   flowEMMi(
#     frame = flowData_transformed_asinh_dn[[1]],
#     ch1 = param[1],
#     ch2 = param[3],
#     sample_size = 1,
#     prior = FALSE,
#     separation = TRUE,
#     max_inits = 10,
#     use_log = FALSE,
#     alpha = .7,
#     img_format = "png",
#     foreground_maxsd = 5000,
#     start_cluster = 2,
#     end_cluster = 12
#   )
```

# Session info


```r
sessionInfo()
```

```
## R version 4.0.1 (2020-06-06)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 17134)
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=Dutch_Belgium.1252  LC_CTYPE=Dutch_Belgium.1252   
## [3] LC_MONETARY=Dutch_Belgium.1252 LC_NUMERIC=C                  
## [5] LC_TIME=Dutch_Belgium.1252    
## 
## attached base packages:
## [1] grid      stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] randomcoloR_1.1.0.1       gtools_3.8.2             
##  [3] mixtools_1.2.0            RcppEigen_0.3.3.7.0      
##  [5] Rcpp_1.0.5                readxl_1.3.1             
##  [7] FlowSOM_1.20.0            igraph_1.2.5             
##  [9] tidyr_1.1.1               RColorBrewer_1.1-2       
## [11] ggcyto_1.16.0             flowWorkspace_4.0.6      
## [13] ncdfFlow_2.34.0           BH_1.72.0-3              
## [15] RcppArmadillo_0.9.900.2.0 cowplot_1.0.0            
## [17] scales_1.1.1              ggplot2_3.3.2            
## [19] dplyr_1.0.1               plyr_1.8.6               
## [21] Phenoflow_1.1.2           foreach_1.5.0            
## [23] flowAI_1.18.5             flowFDA_0.99             
## [25] mclust_5.4.6              multcomp_1.4-13          
## [27] TH.data_1.0-10            MASS_7.3-51.6            
## [29] survival_3.2-3            mvtnorm_1.1-1            
## [31] flowFP_1.46.0             flowViz_1.52.0           
## [33] lattice_0.20-41           flowClean_1.26.0         
## [35] flowCore_2.0.1           
## 
## loaded via a namespace (and not attached):
##   [1] changepoint_2.2.2           ConsensusClusterPlus_1.52.0
##   [3] splines_4.0.1               digest_0.6.25              
##   [5] htmltools_0.5.0             magrittr_1.5               
##   [7] CytoML_2.0.5                cluster_2.1.0              
##   [9] sfsmisc_1.1-7               recipes_0.1.13             
##  [11] Biostrings_2.56.0           gower_0.2.2                
##  [13] RcppParallel_5.0.2          matrixStats_0.56.0         
##  [15] sandwich_2.5-1              cytolib_2.0.3              
##  [17] jpeg_0.1-8.1                colorspace_1.4-1           
##  [19] xfun_0.16                   crayon_1.3.4               
##  [21] jsonlite_1.7.0              hexbin_1.28.1              
##  [23] graph_1.66.0                zoo_1.8-8                  
##  [25] iterators_1.0.12            ape_5.4                    
##  [27] glue_1.4.1                  gtable_0.3.0               
##  [29] ipred_0.9-9                 zlibbioc_1.34.0            
##  [31] XVector_0.28.0              V8_3.2.0                   
##  [33] phyloseq_1.32.0             kernlab_0.9-29             
##  [35] IDPmisc_1.1.20              Rgraphviz_2.32.0           
##  [37] Rhdf5lib_1.10.1             BiocGenerics_0.34.0        
##  [39] bit_4.0.4                   stats4_4.0.1               
##  [41] tsne_0.1-3                  lava_1.6.7                 
##  [43] prodlim_2019.11.13          ellipsis_0.3.1             
##  [45] farver_2.0.3                pkgconfig_2.0.3            
##  [47] XML_3.99-0.5                nnet_7.3-14                
##  [49] caret_6.0-86                labeling_0.3               
##  [51] tidyselect_1.1.0            rlang_0.4.7                
##  [53] reshape2_1.4.4              munsell_0.5.0              
##  [55] cellranger_1.1.0            tools_4.0.1                
##  [57] generics_0.0.2              ade4_1.7-15                
##  [59] evaluate_0.14               biomformat_1.16.0          
##  [61] stringr_1.4.0               yaml_2.2.1                 
##  [63] ModelMetrics_1.2.2.2        knitr_1.29                 
##  [65] purrr_0.3.4                 RBGL_1.64.0                
##  [67] nlme_3.1-148                xml2_1.3.2                 
##  [69] compiler_4.0.1              curl_4.3                   
##  [71] png_0.1-7                   tibble_3.0.3               
##  [73] stringi_1.4.6               Matrix_1.2-18              
##  [75] vegan_2.5-6                 permute_0.9-5              
##  [77] multtest_2.44.0             vctrs_0.3.2                
##  [79] pillar_1.4.6                lifecycle_0.2.0            
##  [81] data.table_1.13.0           R6_2.4.1                   
##  [83] latticeExtra_0.6-29         KernSmooth_2.23-17         
##  [85] gridExtra_2.3               RProtoBufLib_2.0.0         
##  [87] IRanges_2.22.2              codetools_0.2-16           
##  [89] boot_1.3-25                 rhdf5_2.32.2               
##  [91] withr_2.2.0                 S4Vectors_0.26.1           
##  [93] mgcv_1.8-31                 parallel_4.0.1             
##  [95] rpart_4.1-15                timeDate_3043.102          
##  [97] class_7.3-17                rmarkdown_2.3              
##  [99] segmented_1.2-0             Rtsne_0.15                 
## [101] pROC_1.16.2                 Biobase_2.48.0             
## [103] lubridate_1.7.9             base64enc_0.1-3
```
