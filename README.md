# DYNAMO

These are the codes associated to the article Santolini, M., Barabási, A.-L., Recovering perturbation patterns of biological networks from their topology (in preparation).

A test model BIOMD0000000404.xml corresponding to bacterial chemotaxis is provided. 

## Jacobian computation
First compute the Jacobian using compModelSensitivityAndJacobian.m. It is a matlab code and it needs libSBML (http://sbml.org/Software/libSBML/Downloading_libSBML).

## Computation of DYNAMO accuracies
The core of the DYNAMO computation and accuracy evaluation are done in compAccuracies.R (reproduces figure 4b). 

## Robustness to network perturbation
The code corresponding to network perturbation (as in Figure 2g) is compEdgeRemoval.R.

