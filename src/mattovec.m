% Convert a symmetric matrix to a vector 

function v = mattovec(M)
n = size(M, 1);
v = zeros(n*(n+1)/2, 1);
l = 0;
for i = 1:n
   v(l+1:l+i-1) = sqrt(2)*M(1:i-1, i);
   v(l+i) = M(i, i);
   l = l + i;
end
end
