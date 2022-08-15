% Convert a symmetric matrix to a vector 

function v = Mat2Vec(M)
n = size(M, 1);
v = zeros(n*(n+1)/2, 1);
l = 1;
for i = 1:n
   v(l) = M(i, i);
   v(l+1:l+n-i) = sqrt(2)*M(i, i+1:end);
   l = l + n - i + 1;
end
end
