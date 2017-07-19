
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
The source code files and data folder must be in the same folder.


DESCRIPTION OF KEY FILES:
----
- ENN_Demo3.m: This is the main file for the demo. 

- ENNV2Test.m: This is the m function that implements our ENN.V2 method for classification. Given a data set, it outputs the performance of KNN and ENN (ENN.V2 is implemented).  

- buildWeightedKNNMap.m: This is the m function that builds weighted K nearest neighbors distance tree. This function is used in ENNV2Test.m. 


FOR MORE DETAILS, REFER TO OUR PAPER:
B. Tang and H. He, "ENN: Extended Nearest Neighbor Method for Multivariate Pattern Classification," IEEE Computational Intelligence Magazine. 