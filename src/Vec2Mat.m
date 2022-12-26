% Convert a vector to a symmetric matrix

function M = Vec2Mat(v, n)
M = zeros(n, n);
l = 1;
for i = 1:n
   M(i, i) = v(l);
   M(i, i+1:end) = v(l+1:l+n-i)/sqrt(2);
   M(i+1:end, i) = v(l+1:l+n-i)/sqrt(2); 
   l = l + n - i + 1;
end
end