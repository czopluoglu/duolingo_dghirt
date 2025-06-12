# Enhanced Deterministic Gated Hierarchical Item Response Model to Simultaneously Identify Compromised Items and Examinees with Item Preknowledge Using Both Response Accuracy and Response Time Data

**Acknowledgement**

The content of this GitHub repo is a product of a research project funded by Duolingo, Inc. through Competitive Research Grant Program to support topics of interest related to Duolingo's English Test's ongoing research agenda.

**Repository Overview**

For a tutorial-style introduction to the analyses conducted in the paper, please visit:

https://czopluoglu.github.io/duolingo_dghirt

**Directory Structure and Contents**

1. **/data/**:

    - Contains the simulated dataset in both long and wide format.

2. **/script/**:
   
    - R script for simulating the dataset with a mix of honest and fraudulent responses

    - R Script for fitting the model in Stan and extracting model parameters, and evaluating the model performance

    - The Stan model syntax (dghirt.stan) and its compiled version for the cmdstanr package.
  
3. **/docs/**:
   
    - Contains the Rmarkdown file and the rendered HTML for the tutorial page.
