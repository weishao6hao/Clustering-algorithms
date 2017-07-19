clear all;
data=load('can383.txt');
[m n]=size(data);
x_mean=sum(data)/m;
D=pdist2(x_mean,data);
beita=sum((D.^2)')/m;
diata=0.97;
Js=[];
%J=CAA(data,beita,1);
for i=1:5
    m=5*i;
    J=CAA(data,beita,m);
    Js=[Js J];
end
for i=1:5
   p=anova1(Js(:,1:i));
   if 1-p>0.97
       m=i-1;
       break
   end
end
figure(1)
plot(data(:,1),data(:,2),'*');
gamma=m*5;
J=CAA(data,beita,gamma);
[z]=SCA(data,beita,gamma);   %Iterates until it converges

figure(2)
plot3(data(:,1),data(:,2),J,'*');
class=kmeans(z,2);
figure(3)
gscatter(data(:,1),data(:,2),class,'rgbmykc','xd*+pso');