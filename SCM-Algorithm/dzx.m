function output = dzx(X,Z)
[n,d] = size(X);
D = zeros(n);
for i = 1:n
    d2 = zeros(n,1);
    for j = 1:d
        d2 = d2 + (X(i,j) - Z(:,j)).^2;
    end
    D(:,i) = d2;
end
output = D;