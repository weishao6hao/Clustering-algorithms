 clear all;
% X=load('Aggregation.txt');
load data_test.mat
% X=im;
 X=data;
%    X=X(:,1:2);
%    X=normalize(X);
% load data_random5;
%  X=data;
[m n]=size(X);
%  X=double(X);
%  X=normalize(X);
% k=40;J_temp=zeros(m,1);

beta = sum(var(X,1));
%  beta=1;
K = exp(-dzx(X,X)/beta);
%Js=sum(K);
%------------ CCA ------------ 
m = 1;
tol1 = 0.995;
rho = 0;
Js = sum(K.^5);
while rho < tol1
    rm = 5*m;
    rm1 = 5*(m+1);
    Js1 = sum(K.^rm1);
    co = corrcoef(Js,Js1);
    rho = co(2);
    m = m + 1;
    Js = Js1;
    r = rm;
end
Js=Js;
K = round(sqrt(size(X,1)));  
%  K=3;

[rhoArray, distMatrix] = RhoCalculation(X,K);
% [rhoArray, distMatrix] = RhoCalculationAA(data, 2);
rhoArray=Js';
[deltaArray, nneighArray] = DeltaCalculation(distMatrix, rhoArray);
gama = rhoArray.*deltaArray;
% plot(rhoArray, deltaArray,'.')
% plot(sort(gama),'r.')
% hist(gama)
[centerIndexes] = CenterCalculate(X,gama, deltaArray);
[labels] = CluteringBasedOnCenter(rhoArray, nneighArray, centerIndexes);
gscatter(X(:,1),X(:,2),labels,'rgbcmyk','*+<>dox');
figure(2)
plot3(X(:,1),X(:,2),Js,'.');
figure(3)
plot(sort(gama),'.');
