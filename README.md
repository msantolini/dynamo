# DYNAMO

These are the codes associated to the article: 

**Santolini,M. and Barabasi,A.-L. (2018) Predicting perturbation patterns from the topology of biological networks. Proc Natl Acad Sci USA, 169, 201720589.** ([pdf](https://marcsantolini.files.wordpress.com/2018/06/1720589115-full.pdf))

A test model BIOMD0000000404.xml corresponding to bacterial chemotaxis is provided. 

## Jacobian computation
First compute the Jacobian using compModelSensitivityAndJacobian.m. It is a matlab code and it requires libSBML to be run (http://sbml.org/Software/libSBML/Downloading_libSBML).

## Computation of DYNAMO accuracies
The core of the DYNAMO computation and accuracy evaluation are done in compAccuracies.R (reproduces figure 4b). 

## Robustness to network perturbation
The code corresponding to network perturbation (as in Figure 2g) is compEdgeRemoval.R.

