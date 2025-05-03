# R Code and Data to the Article: Class-Focused Variable Importance in Random Forests for Multi-Class Outcomes

Authors: Roman Hornung<sup>1,2,*</sup> and Alexander Hapfelmeier<sup>3</sup>

1. Institute for Medical Information Processing, Biometry and Epidemiology, Faculty of Medicine, LMU Munich, Marchioninistr. 15, Munich, 81377, Germany, ORCID: 0000-0002-6036-1495.
2. Munich Center for Machine Learning (MCML), Munich, Germany.
3. Institute of AI and Informatics in Medicine, TUM School of Medicine and Health, Technical University of Munich, Ismaninger Str. 22, Munich, 81675, Germany, ORCID: 0000-0001-6765-6352.

\* For questions, please contact: hornung@ibe.med.uni-muenchen.de

---

## Program and Platform

- **Program**: R, versions 4.2.2 and 4.4.0.
- The raw results of the simulation study were obtained on a Linux cluster, and the evaluation of the raw results to produce the final results as well as the conduction of the real data analyses was performed on Windows 11.
- Below is the output of the R command `sessionInfo()` on the Linux machine and on the Windows machine. 
  The output specifies which R packages and versions of those packages were used to generate the raw results 
  and to evaluate them.

### sessionInfo() on the Linux cluster

```R
> sessionInfo()
R version 4.2.2 Patched (2022-11-10 r83330)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Debian GNU/Linux 12 (bookworm)

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.21.so

locale:
 [1] LC_CTYPE=de_DE.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=de_DE.UTF-8        LC_COLLATE=de_DE.UTF-8    
 [5] LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=de_DE.UTF-8   
 [7] LC_PAPER=de_DE.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats     graphics  grDevices utils    
[6] datasets  methods   base     

other attached packages:
[1] diversityForest_0.6.0 ranger_0.14.1        
[3] doParallel_1.0.17     iterators_1.0.14     
[5] foreach_1.5.2        

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.10      pillar_1.10.1    compiler_4.2.2  
 [4] ggpubr_0.6.0     tools_4.2.2      lifecycle_1.0.4 
 [7] tibble_3.2.1     gtable_0.3.6     lattice_0.20-45 
[10] pkgconfig_2.0.3  rlang_1.1.1      Matrix_1.5-3    
[13] cli_3.6.1        rstudioapi_0.14  dplyr_1.1.3     
[16] generics_0.1.3   vctrs_0.6.3      grid_4.2.2      
[19] tidyselect_1.2.1 glue_1.8.0       R6_2.6.1        
[22] rstatix_0.7.2    Formula_1.2-5    carData_3.0-5   
[25] ggplot2_3.5.1    purrr_1.0.4      tidyr_1.3.1     
[28] car_3.1-3        magrittr_2.0.3   scales_1.3.0    
[31] backports_1.5.0  codetools_0.2-19 abind_1.4-8     
[34] colorspace_2.1-1 ggsignif_0.6.4   munsell_0.5.1   
[37] broom_1.0.7   
```

### sessionInfo() on the Windows machine

```R
> sessionInfo()
R version 4.4.0 (2024-04-24 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26100)

Matrix products: default


locale:
[1] LC_COLLATE=German_Germany.utf8  LC_CTYPE=German_Germany.utf8   
[3] LC_MONETARY=German_Germany.utf8 LC_NUMERIC=C                   
[5] LC_TIME=German_Germany.utf8    

time zone: Europe/Berlin
tzcode source: internal

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets 
[7] methods   base     

other attached packages:
 [1] OpenML_1.12           cowplot_1.1.3         RColorBrewer_1.1-3   
 [4] tidyr_1.3.1           forcats_1.0.0         stringr_1.5.1        
 [7] xtable_1.8-4          dplyr_1.1.4           gridExtra_2.3        
[10] purrr_1.0.2           patchwork_1.2.0       ggplot2_3.5.1        
[13] ranger_0.17.0         diversityForest_0.6.0

loaded via a namespace (and not attached):
 [1] generics_0.1.3    rstatix_0.7.2     stringi_1.8.4    
 [4] lattice_0.22-6    digest_0.6.35     magrittr_2.0.3   
 [7] evaluate_0.23     BBmisc_1.13       pkgload_1.3.4    
[10] fastmap_1.2.0     jsonlite_1.9.0    Matrix_1.7-0     
[13] backports_1.5.0   Formula_1.2-5     httr_1.4.7       
[16] scales_1.3.0      XML_3.99-0.16.1   abind_1.4-8      
[19] cli_3.6.2         rlang_1.1.3       munsell_0.5.1    
[22] yaml_2.3.8        cachem_1.1.0      withr_3.0.2      
[25] tools_4.4.0       memoise_2.0.1     checkmate_2.3.1  
[28] ggsignif_0.6.4    colorspace_2.1-1  ggpubr_0.6.0     
[31] curl_5.2.1        broom_1.0.7       vctrs_0.6.5      
[34] R6_2.6.1          lifecycle_1.0.4   car_3.1-3        
[37] pkgconfig_2.0.3   pillar_1.10.1     gtable_0.3.6     
[40] glue_1.7.0        data.table_1.15.4 Rcpp_1.0.14      
[43] xfun_0.44         tibble_3.2.1      tidyselect_1.2.1 
[46] knitr_1.46        rstudioapi_0.16.0 htmltools_0.5.8.1
[49] rmarkdown_2.27    carData_3.0-5     compiler_4.4.0
```

---

## General Information and Contents of this Electronic Appendix

### Preliminary Remark
Readers who are not interested in the detailed 
  contents of this electronic appendix, but only in the evaluation of 
  the results or the full reproduction of the simulation study, may skip to 
  the sections "Evaluation of the Results" or "Full Reproduction 
  of the Simulation Study", respectively.

### Contents
- **simulation**: This subfolder contains the R scripts `simulation.R`,
    `simulation_functions.R`, and `simulation_evaluation.R` as well
    as the subfolder `intermediate_results`.
    
    The R script `simulation.R` can be used to run the simulation study 
    using parallel computing on the Linux cluster, producing the 
    raw results (see the section "Full Reproduction of the Simulation Study" 
    below for details).

    The R script `simulation_functions.R` is sourced by `simulation.R` 
    and contains all of the functions required for the simulation study.

    The R script `simulation_evaluation.R` evaluates the raw results of 
    the simulation study to produce all figures and tables associated with 
    the simulation study.

    The subfolder `intermediate_results` contains the intermediate results 
    of the simulation study.
- **application**: This subfolder contains the R scripts `application.R` and `plot_funs.R` as well as the Rda file `gas-drift.Rda`.

    The R script `application.R` performs the real data analyses and produces all corresponding figures and tables.

    The R script `plot_funs.R` contains modified versions of the visualization functions from the diversityForest R package.     These modifications ensure that the resulting figures are suitable for printing in black and white.

    The Rda file `gas-drift.Rda` is the gas-drift dataset downloaded from OpenML (ID: 1476).
- **figures**: This subfolder contains all figures shown in the main paper and in the supplementary material as eps and pdf files, respectively.
- **tables**: This subfolder contains the tables shown in the supplementary material as tex files. Note that when the tables were included in the supplementary material, the tex code of these files was slightly modified for visual reasons (without changing the values in the tables).
    Note further that Table 1 is not included here because this table is related to the design of the simulation study and thus does not show empirical or simulation results.

---

## Evaluation of the Results

The R script `simulation_evaluation.R` contained in the subfolder `simulation` produces all results of the simulation study without the need of re-performing it. This R script reads in a Rda file (stored in the subfolder `simulation/intermediate_results`) that contains the raw results.

The R script `application.R` performs the real data analyses and produces all results of these analyses.

---

## Full Reproduction of the Simulation Study
- As a first step, the folder `vim-mc-rf` this README is contained
  in has to be placed in the home directory (`~/`) of a Linux machine.
  Alternatively, the folder can be placed in a directory of your choice, with the paths 
  in line 3 of the R script `simulation/simulation.R` being changed accordingly.
- An MPI environment is required.
- The R script `simulation/simulation.R` performs the simulation study.

  Note that you will need to change the value of `ncores_used` (i.e., the number of parallel computations) in accordance with your system's resources.

  The script uses the R packages `parallel` and `doParallel`. However, it is also possible to use other parallelization techniques or even sequential computation to reproduce the results. This is because we use a different seed  for each line in the `scenariogrid` data frame created by the `simulation.R` script. Each line in this data frame corresponds to one replication in the simulation study (see `simulation.R` for details). This makes the reproducibility independent of the specific type of parallelization. However, to use a different type of parallelization or to use sequential computation, it is necessary to change the `simulation.R` script accordingly.
