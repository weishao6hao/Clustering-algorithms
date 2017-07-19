close all;
clc;
data=load('can383.txt'); 
[class re]=KMeans(data,5); 
gscatter(data(:,1),data(:,2),class,'rgbcmyk','<>*^.');
