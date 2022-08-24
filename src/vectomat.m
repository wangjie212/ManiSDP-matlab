% Convert a vector to a symmetric matrix

function M = vectomat(v, n)
M = zeros(n, n);
l = 0;
for i = 1:n
   M(1:i-1, i) = v(l+1:l+i-1)/sqrt(2);
   M(i, 1:i-1) = v(l+1:l+i-1)/sqrt(2);
   M(i, i) = v(l+i);
   l = l + i;
end
end