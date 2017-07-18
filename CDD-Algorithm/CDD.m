clear all;
tic;
load data_test3.mat

% data=load('t4.8k.dat');
X=data;
X=normalize(X);
[m n]=size(X);
[D I]=pdist2(X,X,'euclidean','Smallest',m);
I=I';

% Calculate density
k=10;
for i=1:m
    J(i,1)=1/sum(D(2:k+1,i));
end
[r c]=sort(J);
deri=diff(sort(J));

p=find(max(deri(0.05*m:0.95*m,1))==deri);
% p=475;
T=J(c(p,1),1);
figure(1)
plot(X(:,1),X(:,2),'.');
xlabel('X');ylabel('Y')

% title('数据集D4');
figure(2)
plot(deri);
xlabel('样本点');ylabel('密度差分值Δρ')

figure(3)
plot(X(find(J<T),1),X(find(J<T),2),'k>');
hold on;

plot(X(find(J>T),1),X(find(J>T),2),'b*');

legend('噪声点','有用点');

% title('CDD');
xlabel('X');ylabel('Y')
figure(4)
plot(sort(J));
xlabel('样本点');ylabel('密度ρ')
% Remove outliers
class=zeros(m,1);
class(find(J<T))=-1;
lay=1;
ind=[];
rem=find(J>T);
rem1=find(J<T);
data_3N(:,1:2)=[X(rem,1) X(rem,2)];
% Create neighborhood and clustering 
[class nn]=Nature_neighbor(data_3N);
figure(5)
plot(X(rem1,1),X(rem1,2),'rx');
hold on;
gscatter(data_3N(:,1),data_3N(:,2),class,'gbkcymrw','*+>>+*<>p');
xlabel('X');ylabel('Y')
title('CDD(k=10)');

toc