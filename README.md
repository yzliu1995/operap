# An R Package - OPERAP

Here we introduce the OPERAP R package, designed to facilitate the implementation of the lasso tree method and the OPEPA algorithm with pruning. The package offers flexibility in handling both survival outcomes and binary outcomes, allowing for adjustment of non-risk-factor covariates during model fitting. With OPERAP, cancer staging
using ordinal risk factors can be performed, and the staging results can be conveniently saved to a specified file path.

When dealing with time-to-event outcomes, OPERAP also enables the generation of Kaplan-Meier curves, providing a visual assessment of the effectiveness in separating different stages. The package supports various pruning methods, including coarse pruning and fine pruning. Fine pruning can be performed through exhaustive search or quadratic programming constraint.

Additionally, users can define different criteria to determine when to stop pruning. Options include AIC (Akaike Information Criterion), the (integrated) Brier score, the likelihood ratio test, or a predefined total number of stages. OPERAP offers a comprehensive toolkit for efficient and customizable cancer staging analysis.

To illustrate the implementation of cancer staging using our package, we utilize the public METABRIC dataset. The dataset can be found [here](https://www.cbioportal.org/study/summary?id=brca_metabric). It consists of 1,238 patients after removing missing values. Three risk factors, namely Pam50, tumor grade, and neoplasm histologic grade, are employed for cancer staging. Additionally, age at diagnosis can be incorporated as a non-risk-factor covariate. Within our package, the pivotal function for conducting cancer staging is `runOpera()`.

## Getting Started

To install our R package, we will need to use the devtools package. If we haven’t installed `devtools` yet, we can do so by running the following code:

````
install.packages("devtools")
````

Once we have devtools installed, we can proceed with installing our R package:

````
library(devtools)
devtools::install_github("yzliu1995/operap")
````

Make sure to load the package after installation to access the functions and features provided by our package.

````
library(operap)
````

## Data

Once our package is loaded, we can proceed to explore the example dataset included within it. This particular dataset originates from the well-known METABRIC study, which we previously mentioned and discussed.

````
data("bric_bc_os")
````

## Model Fitting

To facilitate cancer staging, there are two initial options available: utilizing OPERA or lasso tree. These methods can be employed to initialize the stages, followed by employing pruning with various stopping rules to refine and reduce the number of stages. In order to provide a comprehensive understanding of this bottom-up approach and how it can be implemented in our R package, we will illustrate the process step by step using our provided example dataset.

### No Pruning

#### OPERA

`runOpera()` is a crucial function to utilize when implementing our R package for cancer staging. Below are explanations for the essential arguments of the function. To perform cancer staging without pruning using OPERA, the following list of arguments needs to be specified:

`ncat`: A vector that indicates the number of levels for each risk factor. In our example dataset, we have three risk factors, with 5 levels in the first one, 4 levels in the second one, and 3 levels in the third one. Therefore, we specify it as `c(5L, 4L, 3L)`.

`dat`: The dataset in the data.frame format. In this case, we specify it as `bric_bc_os`, which contains the example dataset.

`plt`: A flag to determine whether to plot the figures for visualizing the risk categories. By default, it is set to `TRUE`.

`filepath`: The path where all the results and figures should be saved. In our example, it is set as `./results/no pruning/`. If the path does not exist, it will be automatically created. Note that there must be a slash at the end of the file path string.

`riskVariables`: A vector containing the variable names of the risk factors used for staging. In our example, the three variables are `Pam50 + Claudin - low subtype`, `Tumor Stage`, and `Neoplasm Histologic Grade`.

`riskFactors`: A list specifying the levels for each risk factor in ascending order of prognosis. In our example, it is represented as `list(c("LumA", "Normal", "LumB", "Her2", "Basal"), c(1, 2, 3, 4), c(1, 2, 3))`. If a risk factor does not have a total ordering, we need to set the argument `ifE = TRUE` and provide the partial orderings as edges in the network using the `edges` argument.

`riskNames`: A vector containing the names of the risk factors, used to rename them in the resulting figures. In our example, they are `"a - Pam50"`, `"b - Tumor Grade"`, and `"c - Neoplasm Histologic Grade"`.

`TimeN`: The variable name for survival times. In our example, it is `"Overall Survival (Months)"`.

`yN`: The variable name for a binary outcome. In this case, it is NULL since we are using a survival outcome.

`cenN`: The variable name for censoring times. It should be numeric, with zeros indicating censored observations and ones indicating events. In our example, it is `"censoringStatus"`.

`covN`: The variable name(s) for covariates. In our example, it is `"Age at Diagnosis"`.

`withCov`: Specifies whether any covariates need adjustment.

`type`: The type of outcome, either "surv" for survival outcome or "bin" for binary outcome.

`minObs`: The minimum number of patients required in each stage. In our example, any number of patients is allowed.

`legend_size`: The size of the legend in the network of risk categories figure.

`x_pos`: The horizontal position of the legend in the network of risk categories figure.

`yratio`: A ratio proportional to the number of levels in each factor, controlling the vertical position of the legend in the network of risk categories figure.

`yaxis_min`: The minimum value of the y-axis for the Kaplan-Meier curves.

`xpos_bs`: The horizontal position of the label for the Brier score on the x-axis scale.

`x_label`: The label of the x-axis for the Kaplan-Meier curves.

Additionally, there are several other parameters for figure configuration that can be passed to `graphics::title()` and `igraph::plot.igraph()` through our function.

The crucial argument is to set `useOPERA = TRUE` since we specifically intend to employ OPERA as the staging method. Below is an example of the code:

````
r_bric_bc_cp_opera_os <- runOpera(ncat = c(5L, 4L, 3L), 
                                  dat = bric_bc_os, 
                                  plt = T, 
                                  filepath = "./results/no pruning/", 
                                  riskVariables = c("Pam50 + Claudin-low subtype", "Tumor Stage", "Neoplasm Histologic Grade"), 
                                  riskFactors = list(c("LumA", "Normal", "LumB", "Her2", "Basal" ), c(1, 2, 3, 4), c(1, 2, 3)),
                                  riskNames = c("a - Pam50", "b - Tumor Grade", "c - Neoplasm Histologic Grade"), 
                                  useOPERA = T, 
                                  TimeN = "Overall Survival (Months)", 
                                  cenN = "censoringStatus", 
                                  yN = NULL, 
                                  type = "surv", 
                                  covN = "Age at Diagnosis", 
                                  withCov = T, 
                                  minObs = 0,
                                  # parameters for figures
                                  cex.main = 1, 
                                  xpos_bs = 50, 
                                  x_label = "Months",
                                  edge.width = 0.5, 
                                  edge.arrow.size = 0.1, 
                                  vertex.size = 11, 
                                  vertex.frame.color  = NA, 
                                  vertex.label.cex = 0.5, 
                                  vertex.label.color =  "black", 
                                  legend_size = 0.7, 
                                  x_pos = -1.5, 
                                  vertex.color = "lightblue")
````


If users prefer each stage to have a minimum number of patients (e.g., 30), they can modify the argument `minObs` as `minObs = 30`. This adjustment enables coarse pruning to select the optimal staging result, and ensures that each stage consists of at least `30` patients.

In the example above, we assume that each risk factor has a total ordering. However, our package can also deal with some risk factor(s) with only partial ordering. The key is to set the correct edges representing the partial ordering.

For example, if `Pam50` is a partially ordered risk factor with the partial ordering as `LumB <= Her2 <= Basal`, and no ordering for either `LumA` or `Normal` with respect to other levels. The network of all three risk factors will be comprised of one sub-network associated with levels from `LumB <= Her2 <= Basal`, and the other two sub-networks associated with `LumA` and `Normal` respectively. Below is the code that shows how to define the edges in this scenario. To perform the cancer staging, users only need to add `ifE = TRUE` and `edges = edges`. Note that there is no need to change the value passed to `riskFactors` as we still want `a1` to represent `LumA`, `a2` to represent `Normal`, and so on.

````
# One sub-network associated with levels from `LumB <= Her2 <= Basal`
sub_1 <- edgesHasse(ncat = c(3L, 4L, 3L), e =  c())
for(i in 1:3){
  sub_1 <- gsub(paste0("a", i), paste0("e", 2+i), sub_1)
} 
sub_1 <- gsub("e", "a", sub_1)

# The other two sub-networks associated with `LumA` and `Normal`
sub_2 <- edgesHasse(ncat = c(1L, 4L, 3L), e =  c())
sub_3 <-gsub("a1", "a2", sub_2)

edges <- c(sub_2, sub_3, sub_1)
````


#### Lasso Tree

The crucial argument to note is to set `useLassoT = TRUE` instead of `useOPERA = TRUE`, as we have chosen to utilize the lasso tree as the staging method. With this change, all other steps and procedures can be followed in a similar manner to when we used OPERA. Moreover, the tolerance accuracy for convergence can be adjusted by setting the value of `eps_lasso`, which controls the speed at which the algorithm converges. Additionally, during parameter tuning, we have the option to use the AIC instead of the BIC by specifying `useBIC = FALSE`.

### Pruning

To enable pruning, we need to specify the argument `usePruning = TRUE`. There are four different stopping rules available: the likelihood ratio test with a pre-specified Type I error rate $\alpha$, the AIC, the (integrated) Brier score, or a predefined total number of stages. In our package, if we choose the likelihood ratio test, we can set the argument `useLRT = TRUE`, and the corresponding $\alpha$ (e.g., 0.01) can be specified using `threshold = 0.01`. Using the AIC, the (integrated) Brier score, or a predefined total number of stages (e.g., 5) can be achieved by setting `useAIC = TRUE`, `useIbs = TRUE`, or `prefix_stage = 5`, respectively.

We recommend using the likelihood ratio test with an $\alpha = 0.01$ when the final number of stages is unknown, as it has shown superior performance in simulation studies. Alternatively, if the number of stages is known in advance, using a predefined number is a suitable option. There are three fundamental pruning methods available: coarse pruning, fine pruning using exhaustive search, and fine pruning using quadratic constraint. In our package, coarse pruning can be implemented by setting `coarse_pruning = TRUE`, fine pruning using exhaustive search by setting `fine_pruning = TRUE`, and fine pruning using quadratic constraint by setting `fine_pruning_quad = TRUE`. It is important to choose only one of the three fundamental pruning methods, depending on the available computational resources.

We recommend fine pruning using exhaustive search due to its superior performance in simulation studies. However, it can be computationally expensive when the total number of risk categories is high. Fine pruning using quadratic constraint and coarse pruning offer decent accuracy comparable to fine pruning using exhaustive search but with less computational burden. Therefore, we suggest choosing either fine pruning using exhaustive search, fine pruning using quadratic constraint, or coarse pruning based on the computational feasibility in the specific scenario.

To demonstrate the implementation of pruning using our package, we will use the likelihood ratio test with an $\alpha = 0.01$ and apply fine pruning using exhaustive search. The same steps can be followed for other methods by specifying the corresponding arguments, as discussed earlier. The code snippet is provided below:

````
r_bric_bc_fp_opera_os <- runOpera(ncat = c(5L, 4L, 3L), 
                                  dat = bric_bc_os, 
                                  plt = T, 
                                  filepath = "./results/pruning/LRT/", 
                                  riskVariables = c("Pam50 + Claudin-low subtype", "Tumor Stage", "Neoplasm Histologic Grade"), 
                                  riskFactors = list(c("LumA", "Normal", "LumB", "Her2", "Basal" ), c(1, 2, 3, 4), c(1, 2, 3)), 
                                  riskNames = c("a - Pam50", "b - Tumor Grade", "c - Neoplasm Histologic Grade"),
                                  useOPERA = T, 
                                  # adds the arguments for pruning
                                  usePruning = T,
                                  fine_pruning = T,
                                  useLRT = T,
                                  threshold = 0.01,
                                  # Other arguments
                                  TimeN = "Overall Survival (Months)", 
                                  cenN = "censoringStatus", 
                                  yN = NULL, 
                                  type = "surv", 
                                  covN = "Age at Diagnosis", 
                                  withCov = T, 
                                  minObs = 0,
                                  # parameters for figures
                                  cex.main = 1, 
                                  xpos_bs = 50, 
                                  x_label = "Months",
                                  edge.width = 0.5, 
                                  edge.arrow.size = 0.1, 
                                  vertex.size = 11, 
                                  vertex.frame.color  = NA, 
                                  vertex.label.cex = 0.5, 
                                  vertex.label.color =  "black", 
                                  legend_size = 0.7, 
                                  x_pos = -1.5, 
                                  vertex.color = "lightblue")
````

## Acknowledgments/References/Citations

* [opera](https://github.com/WangTJ/opera)
* [glidars](https://github.com/WangTJ/glidars)
* [README.md template](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2#file-readme-template-md) 
