#' A \link[SingleCellExperiment]{SingleCellExperiment} object
#'
#' This dataset contains all the information used for the TGFB case studies
#'
#' #' @references
#' Jos\'{e}  L.  McFaline-Figueroa,  Andrew  J.  Hill,  Xiaojie  Qiu,
#' Dana  Jackson,  Jay  Shendure,  and  Cole Trapnell.
#' *A pooled single-cell genetic screen identifies regulatory checkpoints in *
#' *the continuum of the epithelial-to-mesenchymal transition.*
#' Nature Genetics, 51(9):1389–1398, sep 2019.  ISSN 15461718.
#' doi:  10.1038/s41588-019-0489-5.
#' @source
#' See the TGFB vignette to generate this data
#' @usage data(tgfb)
"tgfb"

#' A \link[data.frame]{data.frame} object
#'
#' This contains the output of running the \link[tradeSeq]{conditionTest} on
#' the \link{tgfb} dataset after running \link{monocle3}.
#'
#' #' @references
#' Jos\'{e}  L.  McFaline-Figueroa,  Andrew  J.  Hill,  Xiaojie  Qiu,
#' Dana  Jackson,  Jay  Shendure,  and  Cole Trapnell.
#' *A pooled single-cell genetic screen identifies regulatory checkpoints in *
#' *the continuum of the epithelial-to-mesenchymal transition.*
#' Nature Genetics, 51(9):1389–1398, sep 2019.  ISSN 15461718.
#' doi:  10.1038/s41588-019-0489-5.
#' @source
#' See the TGFB vignette to generate this data
#' @usage data(condRes_monocle)
"condRes_monocle"

#' A \link[SingleCellExperiment]{SingleCellExperiment} object
#'
#' This dataset contains all the information used for the TCDD case studies
#'
#' @references
#' Rance Nault,  Kelly A. Fader,  Sudin Bhattacharya,  and Tim R. Zacharewski.
#' *Single-Nuclei RNA Sequencing Assessment of the Hepatic Effects of*
#' * 2,3,7,8-Tetrachlorodibenzo-p-dioxin.CMGH*,
#' 11(1):147–159, jan 2021.  ISSN 2352345X.
#' doi:  10.1016/j.jcmgh.2020.07.012.
#' @source
#' See the TCDD vignette to generate this data
#' @usage data(tcdd)
"tcdd"

#' A \link[SingleCellExperiment]{SingleCellExperiment} object
#'
#' This dataset contains all the information used for the KRAS case studies
#'
#' @references
#' Jenny Y. Xue, Yulei Zhao, Jordan Aronowitz, Trang T. Mai, Alberto Vides,
#' Besnik Qeriqi, Dongsung Kim, Chuanchuan Li, Elisa de Stanchina, Linas Mazutis,
#'  Davide Risso, and Piro Lito.
#' *Rapid non-uniform adaptation to conformation-specific KRAS(G12C) inhibition.*
#' Nature, 577(7790):421–425,jan  2020. ISSN  14764687. doi: 10.1038/s41586-019-1884-x
#' @source
#' See the KRAS vignette to generate this data
#' @usage data(kras)
"kras"

#' A \link[SingleCellExperiment]{SingleCellExperiment} object
#'
#' This dataset contains all the information used for the fibrosis case study
#'
#' @references
#' Arun  C.  Habermann,  Austin  J.  Gutierrez,  Linh  T.  Bui,  
#' Stephanie  L. Yahn,  Nichelle  I.  Winters, Carla L. Calvi, Lance Peter, 
#' Mei I. Chung, Chase J. Taylor, Christopher Jetter, Latha Raju, Jamie Roberson,
#'  Guixiao Ding, Lori Wood, Jennifer M.S. Sucre, Bradley W. Richmond, Ana P. Serezani,
#'  Wyatt J. McDonnell,  Simon B. Mallal,  Matthew J. Bacchetta,  James E. Loyd,
#'  Ciara M. Shaver,Lorraine B. Ware,  Ross Bremner,  Rajat Walia,  
#'  Timothy S. Blackwell,  Nicholas E. Banovich,  and Jonathan A. Kropski.
#' *Single-cell rna sequencing reveals profibrotic roles of distinct epithelial and mesenchymal  lineages  in  pulmonary  fibrosis*
#' Science  Advances,  6,  2020.   ISSN  23752548.   doi:10.1126/sciadv.aba1972.
#' @source
#' See the fibrosis vignette to generate this data
#' @usage data(fibrosis)
"fibrosis"

#' Results on the one bifurcation and two conditions datasets
#'
#' This dataset contains all the results from the simulation on the datasets
#' with one bifurcation and two conditions
#'
#' @source
#' Output of the simulations
#' @usage data(fork)
#' @format Results on 2000 datasets
"fork"

#' Results on the two bifurcations and two conditions datasets
#'
#' This dataset contains all the results from the simulation on the datasets
#' with two bifurcations and two conditions
#'
#' @source
#' Output of the simulations
#' @usage data(tree)
#' @format Results on 2000 datasets
"tree"

#' Results on the one bifurcation and three conditions datasets
#'
#' This dataset contains all the results from the simulation on the datasets
#' with one bifurcation and three conditions
#'
#' @source
#' Output of the simulations
#' @usage data(complex)
#' @format Results on 900 datasets
"complex"

#' Results on the five lineages datasets
#'
#' This dataset contains all the results from the simulation on the datasets
#' with five lineages
#'
#' @source
#' Output of the simulations
#' @usage data(five_lins)
#' @format Results on 900 datasets
"five_lins"

#' A \link[SingleCellExperiment]{SingleCellExperiment} example of a one bifurcation and two conditions dataset
#'
#' Simulated dataset with one bifurcation and two conditions
#'
#' @source
#' Created by running "fork_sce <- create_bifurcating_simu(multiplier = .5)"
#' @usage data(fork_sce)
#' @format A \link[SingleCellExperiment]{SingleCellExperiment} object
"fork_sce"

#' A \link[SingleCellExperiment]{SingleCellExperiment} example of a two bifurcations and two conditions dataset
#'
#' Simulated dataset with two bifurcations and two conditions
#'
#' @source
#' Created by running "tree_sce <- create_consecutive_bifurcating_simu(multiplier = .5)"
#' @usage data(tree_sce)
#' @format A \link[SingleCellExperiment]{SingleCellExperiment} object
"tree_sce"

#' A \link[SingleCellExperiment]{SingleCellExperiment} example of a one bifurcation and three conditions dataset
#'
#' Simulated dataset with one bifurcation and three conditions
#'
#' @source
#' Created by running "complex_sce <- create_bifurcating_three_conditions_simu(multiplier = .5)"
#' @usage data(complex_sce)
#' @format A \link[SingleCellExperiment]{SingleCellExperiment} object
"complex_sce"

#' A \link[SingleCellExperiment]{SingleCellExperiment} example of 5 lineages and
#'  two conditions dataset
#'
#' Simulated dataset with 5 lineages and two conditions
#'
#' @source
#' Created by running "fivelin_sce <- create_5_lineages_simu(multiplier = .5)"
#' @usage data(fivelin_sce)
#' @format A \link[SingleCellExperiment]{SingleCellExperiment} object
"fivelin_sce"
