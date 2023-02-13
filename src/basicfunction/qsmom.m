% Generate the second-order moment-SDP relaxation for the quartic sphere program:
% Min  coe'*[x]_4
% s.t. x^2 = 1.
% Output sedumi format data.

function [At, b, c, mb] = qsmom(n, coe)
basis = get_basis(n, 2);
mb = size(basis, 2);
sp = get_basis(n, 4);
lsp = size(sp, 2);
mm = cell(lsp,1);
for i = 1:mb
    for j = i:mb
         bi = basis(:,i) + basis(:,j);
         ind = bfind(sp, lsp, bi, n);
         mm{ind} = [mm{ind} [i;j]];
    end
end
ncons = mb*(mb+1)/2 - lsp + mb + 1;
row = [1];
col = [1];
val = [1];
b = sparse(ncons, 1);
b(1) = 1;
l = 2;
loa = cell(lsp,1);
for i = 1:lsp
     loa{i} = zeros(2*size(mm{i},2),1);
     for j = 1:size(mm{i},2)
          loa{i}(2*j-1:2*j) = [(mm{i}(2,j)-1)*mb + mm{i}(1,j);(mm{i}(1,j)-1)*mb + mm{i}(2,j)];
     end
end
for i = 1:mb
    for k = 1:n
        s1 = 0; 
        temp = zeros(n,1);
        temp(k) = 2;
        ind1 = bfind(sp, lsp, basis(:,i)+temp, n);
        for j = 1:size(mm{ind1},2)
            if mm{ind1}(1,j) == mm{ind1}(2,j)
               row = [row; loa{ind1}(2*j)];
               s1 = s1 + 1;
            else
               row = [row; loa{ind1}(2*j-1:2*j)];
               s1 = s1 + 2;
            end
        end
        col = [col; l*ones(s1,1)];
        val = [val; 1/s1*ones(s1,1)];
    end
    ind2 = bfind(sp, lsp, basis(:,i), n);
    s2 = 0;
    for j = 1:size(mm{ind2},2)
        if mm{ind2}(1,j) == mm{ind2}(2,j)
           row = [row; loa{ind2}(2*j)];
           s2 = s2 + 1;
        else
           row = [row; loa{ind2}(2*j-1:2*j)];
           s2 = s2 + 2;
         end
    end
    col = [col; l*ones(s2,1)];
    val = [val; -1/s2*ones(s2,1)];
    l = l + 1;
end

for i = 1:lsp
     [~, idx] = max(mm{i}(1,:));
     for j = 1:size(mm{i},2)
         if j ~= idx
             if mm{i}(1,idx) == mm{i}(2,idx)
                 row = [row; loa{i}(2*idx)];
                 col = [col; l];
                 val = [val; 1];
             else
                 row = [row; loa{i}(2*idx-1:2*idx)];
                 col = [col; l; l];
                 val = [val; 0.5; 0.5];
             end
             if mm{i}(1,j) == mm{i}(2,j)
                 row = [row; loa{i}(2*j)];
                 col = [col; l];
                 val = [val; -1];
             else
                 row = [row; loa{i}(2*j-1:2*j)];
                 col = [col; l; l];
                 val = [val; -0.5; -0.5];
             end
             l = l + 1;
         end
     end
end
At = sparse(row,col,val,mb^2,ncons);

row = [];
col = [];
val = [];
for i = 1:lsp
   s = 0;
   for k = 1:size(mm{i},2)
       if mm{i}(1,k) == mm{i}(2,k)
          row = [row; mm{i}(1,k)];
          col = [col; mm{i}(2,k)];
          s = s + 1;
       else
          row = [row; mm{i}(1,k); mm{i}(2,k)];
          col = [col; mm{i}(2,k); mm{i}(1,k)];
          s = s + 2;
       end
   end
   val = [val; coe(i)/s*ones(s,1)];
end    
C = sparse(row,col,val,mb,mb);
c = C(:);
end