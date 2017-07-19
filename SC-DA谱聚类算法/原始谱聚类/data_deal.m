%后续数据处理
clear all
clc
load IDX_all
load information_dataset
number_loop=20;
rate_correct=zeros(number_loop,1);%记录正确率
CE=zeros(number_loop,1);%记录错误率
NMI=zeros(number_loop,1);%记录NMI值
RI=zeros(number_loop,1);%记录RI值
FI=zeros(number_loop,1);%记录FI值
ARI=zeros(number_loop,1);%记录ARI值
KaI=zeros(number_loop,1);%记录KaI值
for i=1:number_loop
    [corrate]=CorrectRate(X_all,IDX_all(:,i),number_sample,k);%计算正确率及原始中心和新的中心的位置的对应
    value_ce=1-corrate;%计算错误率
    [value_nmi]=NormalizedMutualInformation(X_all,IDX_all(:,i),number_sample,k);%计算聚类结果的NMI值
    [value_ri]=RandIndexM(X_all,IDX_all(:,i),number_sample,k);%计算聚类结果的RandIndex值
    [value_fi]=FIndexM(X_all,IDX_all(:,i),number_sample,k);%计算聚类结果的FIndex值
    [value_ari]=AdjustedRandIndexM(X_all,IDX_all(:,i),number_sample,k);%计算聚类结果的AdjustedRandIndex值
    [value_kai]=KappaIndex(X_all,IDX_all(:,i),number_sample,k);%计算聚类结果的KappaIndex值
    rate_correct(i,1)=corrate;
    CE(i,1)=value_ce;
    NMI(i,1)=value_nmi;
    RI(i,1)=value_ri;
    FI(i,1)=value_fi;
    ARI(i,1)=value_ari;
    KaI(i,1)=value_kai;
end
%均值+方差
mean_rate_correct=mean(rate_correct);%计算平均正确率
var_rate_correct=var(rate_correct);
mean_CE=mean(CE);%计算平均错误率
var_CE=var(CE);
mean_NMI=mean(NMI);%计算平均NMI值
var_NMI=var(NMI);
mean_RI=mean(RI);%计算平均RI值
var_RI=var(RI);
mean_FI=mean(FI);%计算平均FI值
var_FI=var(FI);
mean_ARI=mean(ARI);%计算平均ARI值
var_ARI=var(ARI);
mean_KaI=mean(KaI);%计算平均KaI值
var_KaI=var(KaI);