function [centerIndexes] = OutwardStatTestCenterDetection(gama, delta)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
[data,sInds] = sort(gama, 'descend');
N = length(data);
X = data;
R = data(1:(N-1))./data(2:N);
% R
% bins = hist(X);
% for i = 1:length(bins)-1
%     if bins(i) > bins(i+1)
%         break;
%     end
% end
% 
% T = N - sum(bins([1:i]));
% K = round(0.05*T);
% plot(R,'.')
% if N > 800
%     T = round(0.9*N);
%     K = round(0.1*T);
% %     T = round(0.8*N);
% %     K = round(0.05*T);    
% else
%     T = round(0.9*N);
%     K = 150;
% end

%==================
if N > 800
%     T = round(0.9*N);
%     K = round(0.05*T);
    T = round(0.95*N);
    K = round(0.1*T);    
else
    T = round(0.5*N);
    K = round(0.02*T);
end
%==================
% K = 150;
% K = 200
% alphaK = 1/((1/T)*sum(log(X([1:T]))) - log(X(T+1)));

% K = 150;
% bins = hist(X,10);
% for i = 1:length(bins)-1
%     if bins(i) > bins(i+1)
%         break;
%     end
% end
% 
% T = N - sum(bins([1:i]));
% K = 20;
% T = round(0.5*N);
%     K = 0.5*T; 

% bins = hist(X);
% T = N - sum(bins([1:3]));

% T = round(0.6*N);
% % T = 260;
% K = round(0.3*T);
alphaK = 1/((K/(T - K +1))*log(X(K+1)) - (T/(T-K + 1))*log(X(T+1)) + (1/(T - K + 1))*sum(log(X([(K+1):T]))));
% r = (1 - (1 -delta)^(1/K))^
for i = K:-1:1
    ri = (1 - (1 - delta)^(1/K))^(-1/(alphaK*i));
%     ri
    if R(i) >= ri
        break;
    end
end
 i
%  alphaK
% i
 if i < 2
     i = 2;
 end
centerIndexes = sInds([1:i]);
end

