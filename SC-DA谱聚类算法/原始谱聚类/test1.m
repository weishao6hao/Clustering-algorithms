%实现NJW算法
clear all
clc
number_loop=1;%循环的次数
[X,number_sample,k,class_original]=InputData('D:\参考论文\STClustering Algorithm\flame.txt');%输入数据，得到样本个数、类别数
% k=ceil(k);
X_all=[X;class_original];%矩阵X_all表示数据本身的特征值+类别标号特征
Standard_Deviation=0.1;%人工赋值标准差
rate_correct=zeros(number_loop,1);%记录正确率
number_correct=zeros(number_loop,k);%记录正确样本数
time_run=zeros(number_loop,1);%记录运行时间
CE=zeros(number_loop,1);%记录错误率
NMI=zeros(number_loop,1);%记录NMI值
RI=zeros(number_loop,1);%记录RI值
FI=zeros(number_loop,1);%记录FI值
ARI=zeros(number_loop,1);%记录ARI值
KaI=zeros(number_loop,1);%记录KaI值
IDX_all=zeros(number_sample,number_loop);%记录每次的聚类结果
for i=1:number_loop  
    i
    t1=clock;
    [A]=AffinityMatrix(X,number_sample,k,Standard_Deviation);%step1：根据数据集X，得到相似矩阵A
    [D,L]=LaplacianMatrix(A,number_sample);%step2：定义度矩阵D和利用公式得到矩阵L
    [M]=EigenvectorMatrix(L,number_sample,k);%step3：找到L矩阵的k个最大特征向量，形成M集合，且保证向量彼此正交
    [Y]=Normaliz(M,number_sample,k);%step4：对集合M进行归一化，得到矩阵Y
    [center,index_class,distance_sample,number_sample_class]=K_Means(Y,number_sample,k);%step5：对Y利用kmeans算法进行聚类，得到聚类结果;step6：对Y的分类结果直接映射到原数据集X上
    t2=clock;
    time_run(i)=etime(t2,t1);
    [corrate,numcorrect,ON,number_class]=CorrectRate(X_all,index_class,number_sample,k);%计算正确率及原始中心和新的中心的位置的对应
    value_ce=1-corrate;%计算错误率
    [value_nmi]=NormalizedMutualInformation(X_all,index_class,number_sample,k);%计算聚类结果的NMI值
    [value_ri]=RandIndexM(X_all,index_class,number_sample,k);%计算聚类结果的RandIndex值
    [value_fi]=FIndexM(X_all,index_class,number_sample,k);%计算聚类结果的FIndex值
    [value_ari]=AdjustedRandIndexM(X_all,index_class,number_sample,k);%计算聚类结果的AdjustedRandIndex值
    [value_kai]=KappaIndex(X_all,index_class,number_sample,k);%计算聚类结果的KappaIndex值  
    rate_correct(i,1)=corrate;
    number_correct(i,:)=numcorrect;
    CE(i,1)=value_ce;
    NMI(i,1)=value_nmi;
    RI(i,1)=value_ri;
    FI(i,1)=value_fi;
    ARI(i,1)=value_ari;
    KaI(i,1)=value_kai;
    IDX_all(:,i)=index_class;
end
%多次实验结果的处理及结果的评估
%ShowFigure(X,Y,class_original,index_class,number_sample,k);%分类结果的图像演示
%整理数据
mean_time_run=sum(time_run)/number_loop;%计算平均运行时间
mean_rate_correct=sum(rate_correct)/number_loop;%计算平均正确率
meaan_CE=sum(CE)/number_loop;%计算平均错误率
mean_NMI=sum(NMI)/number_loop;%计算平均NMI值
mean_RI=sum(RI)/number_loop;%计算平均RI值
mean_FI=sum(FI)/number_loop;%计算平均FI值
mean_ARI=sum(ARI)/number_loop;%计算平均ARI值
mean_KaI=sum(KaI)/number_loop;%计算平均KaI值

 gscatter(X(1,:),X(2,:),index_class,'rgbcmyk','<>^*.')
%存储聚类结果
save IDX_all IDX_all;
save information_dataset X_all number_sample k