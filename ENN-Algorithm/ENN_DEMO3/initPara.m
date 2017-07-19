ENN.mode = 'MK';   %% The mode of Monte-Carlo which randomly choose some
%percent of training data and testing data to test the classification performance.
% ENN.mode = 'Standard';  %% The mode of Standard which use all the
% training data and testing data from the Data where the ENN.separateIndex
% has to be existed to separtate the training data and testing data.

ENN.frFlag = false; 
ENN.frMode = 'PCA';  %% The method of PCA to reduce features.

switch datasetName
    case 'haberman'
        ENN.normalizeFlag = false; 
    case 'vowel'
        ENN.normalizeFlag = false; 
    case 'wine'
        ENN.normalizeFlag = true;   
    case 'breast-cancer'
        ENN.normalizeFlag = false;
    case 'breasttissue'
        ENN.normalizeFlag = true;
    case 'ILPD'
        ENN.normalizeFlag = true;
    case 'pima-indians-diabetes'
        ENN.normalizeFlag = true;
    case 'knowledge'
        ENN.normalizeFlag = false;
    case 'Vertebral'
        ENN.normalizeFlag = true;
    case 'magic'
        ENN.normalizeFlag = true;
    case 'mnistData'
        ENN.normalizeFlag = false;
        ENN.mode = 'Standard';
    case 'spamdata'
        ENN.normalizeFlag = true;
    otherwise 
        disp('error: no match title!');
        return
end







