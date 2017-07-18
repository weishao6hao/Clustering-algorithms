clear all;
% data=load('t8.8k.dat');
load data_test3
data=normalize(data);
k=6;
tic

[class,type]=dbscan(data,k,0.2);
gscatter(data(:,1),data(:,2),class,'rbgmykcg','xd*+pso');
xlabel('X');
ylabel('Y');
title('DBSCAN(MinPts=6,Eps=0.2)');
toc

