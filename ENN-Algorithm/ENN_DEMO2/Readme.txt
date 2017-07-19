
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Readme.txt
%% Date: 11/08/2014
%% By Bo Tang (btang@ele.uri.edu) and Haibo He (he@ele.uri.edu)
%% For any questions and/or comments for this code/paper, please feel free to contact Prof. Haibo He, Electrical Engineering, University of Rhode Island, 
%% Email: he@ele.uri.edu
%% Web:  http://www.ele.uri.edu/faculty/he/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DESCRIPTION
----

This code package is the implementation of the proposed Extended Nearest Neighbor (ENN) method for pattern classification. 


DEPENDENCIES
----
- MATLAB 


TO RUN
----
The source code files and data folder must be in the same folder. To run the demos, you have to choose the dataset (the parameter of "datasetName" in the m file).


DESCRIPTION OF KEY FILES:
----
- ENN_Demo2.m: This is the main file for the demo. 

- ENNTest.m: This is the m function that implements our ENN method for classification. Given the training data and test data, it outputs the prediction performance of both KNN and ENN rule.  

- buildWeightedKNNMap.m: This is the m function that builds weighted K nearest neighbors distance tree. This function is used in ENNTest.m. 


FOR MORE DETAILS, REFER TO OUR PAPER:
B. Tang and H. He, “ENN: Extended Nearest Neighbor Method for Multivariate Pattern Classification,” IEEE Computational Intelligence Magazine. 