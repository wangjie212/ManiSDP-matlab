function basis = get_basis(n, d)
lb = nchoosek(n+d, d);
basis = zeros(n, lb);
i = 0;
t = 1;
while i < d + 1
    t = t + 1;
    if basis(n, t-1) == i
       if i < d
          basis(1, t) = i + 1;
       end
       i = i + 1;
    else
        j = 1;
        while basis(j, t-1) == 0
            j = j + 1;
        end
        basis(:, t) = basis(:, t-1);
        if j == 1 
            basis(1, t) = basis(1, t) - 1;
            basis(2, t) = basis(2, t) + 1;
        else
            basis(1, t) = basis(j, t) - 1;
            basis(j, t) = 0;
            basis(j+1, t) = basis(j+1, t) + 1;
         end
    end
end
end