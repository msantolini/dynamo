# DYNAMO

These are codes associated to the article: 

**Santolini,M. and Barabasi,A.-L. (2018) Predicting perturbation patterns from the topology of biological networks. Proc Natl Acad Sci USA, 169, 201720589.** ([pdf](https://marcsantolini.files.wordpress.com/2018/06/1720589115-full.pdf))

A test model BIOMD0000000404.xml corresponding to bacterial chemotaxis is provided. 

## Jacobian computation
First compute the Jacobian using compModelSensitivityAndJacobian.m. It is a Matlab code and it requires the library libSBML to be installed (http://sbml.org/Software/libSBML/Downloading_libSBML#MATLAB). In the paper, SBML Toolbox 4.1.0 was used. The key libSBML functions that are used in our context are in the subfolder codes_libSBML.

## Computation of DYNAMO accuracies
The core of the DYNAMO computation and accuracy evaluation are done in compAccuracies.R (reproduces figure 4b). 

## Comparison between model accuracy and network features
The code plotAccuracyVsFeatures.R computes the association between model accuracy and network features once DYNAMO has been used for several BioModels (reproduces figure 2G). 



