% Generate the second-order moment-SDP relaxation for the sparse binary quadratic program:
% inf  c1'*[x_{I_1}]_2 + ... + ct'*[x_{I_t}]_2
% s.t. x_i^2 = 1, i = 1,...,n.
% Output sedumi format data.

function [At, b, c, K] = bqpmom_sparse(n, I, coe)
t = length(I);
mb = zeros(t, 1);
sp = [];
for i = 1:t
    basis{i} = get_basis(n, 2, I{i});
    ind = true(size(basis{i}, 2), 1);
    for j = 1:size(basis{i}, 2)
        if sum(basis{i}(:,j) > 1) > 0
            ind(j) = false;
        end
    end
    basis{i} = basis{i}(:,ind);
    mb(i) = size(basis{i}, 2);
    temp = get_basis(n, 4, I{i});
    ind = true(size(temp, 2), 1);
    for j = 1:size(temp, 2)
        if sum(temp(:,j) > 2) > 0 || sum(mod(temp(:,j),2)) == 0
           ind(j) = false;
        end
    end
    sp = [sp temp(:,ind)];
end
sp = unique(sp', 'rows');
sp = sortrows(sp)';
lsp = size(sp, 2);
mm = cell(lsp,1);
for k = 1:t
    for i = 1:mb(k)
        for j = i+1:mb(k)
            bi = basis{k}(:,i) + basis{k}(:,j);
            ind = nbfind(sp, lsp, bi, n);
            mm{ind} = [mm{ind} [i;j;k]];
        end
    end
end
mc = zeros(t, 1);
for k = 1:t
    mc(k) = length(I{k});
end
ncons = sum(mb.*(mb+1)/2) - lsp + sum(mc.*(mb-1)) - sum(mb) + t;
row = [1];
col = [1];
val = [1];
b = sparse(ncons, 1);
b(1) = 1;
l = 2;
for k = 1:t
    ss = sum(mb(1:k-1).^2);
    if k == 1
        init = 2;
    else
        init = 1;
    end
    for i = init:mc(k)+1
        row = [row; 1; ss+(i-1)*mb(k)+i];
        col = [col; l; l];
        val = [val; 0.5; -0.5];
        l = l + 1;
    end
end
ind = [1:n];
for k = 1:t
    ss = sum(mb(1:k-1).^2);
    for i = mc(k)+2:mb(k)
        temp = ind(basis{k}(:,i)==1);
        cc = [bfind(I{k}, mc(k), temp(1), 1)+1; bfind(I{k}, mc(k), temp(2), 1)+1];
        row = [row; ss+(cc(1)-1)*mb(k)+cc(1); ss+(i-1)*mb(k)+i; ss+(cc(2)-1)*mb(k)+cc(2); ss+(i-1)*mb(k)+i];
        col = [col; l; l; l+1; l+1];
        val = [val; 0.5; -0.5; 0.5; -0.5];
        l = l + 2;
    end
end
loa = cell(lsp,1);
for i = 1:lsp
     loa{i} = zeros(2*size(mm{i},2),1);
     for j = 1:size(mm{i},2)
          temp = sum(mb(1:mm{i}(3,j)-1).^2);
          loa{i}(2*j-1:2*j) = [temp+(mm{i}(2,j)-1)*mb(mm{i}(3,j))+mm{i}(1,j);temp+(mm{i}(1,j)-1)*mb(mm{i}(3,j))+mm{i}(2,j)];
     end
end
for q = 1:t
    for k = 1:mc(q) 
        for i = 2:mb(q)
            if basis{q}(I{q}(k), i) == 0
                bi = basis{q}(:,i);
                bi(I{q}(k)) = 2;
                ind1 = nbfind(sp, lsp, bi, n);
                ind2 = nbfind(sp, lsp, basis{q}(:,i), n);
                row = [row; loa{ind1}; loa{ind2}];
                col = [col; l*ones(length(loa{ind1})+length(loa{ind2}),1)];
                if length(loa{ind1}) < length(loa{ind2})
                    val = [val; ones(length(loa{ind1}),1); -length(loa{ind1})/length(loa{ind2})*ones(length(loa{ind2}),1)];
                else
                    val = [val; length(loa{ind2})/length(loa{ind1})*ones(length(loa{ind1}),1); -ones(length(loa{ind2}),1)];
                end
                l = l + 1;          
            end
        end
    end
end

for i = 1:lsp
     [~, idx] = max(mm{i}(1,:));
     for j = 1:size(mm{i},2)
         if j ~= idx
             row = [row; loa{i}(2*idx-1:2*idx); loa{i}(2*j-1:2*j)];
             col = [col; l; l; l; l];
             val = [val; 0.5; 0.5; -0.5; -0.5];
             l = l + 1;
         end
     end
end
At = sparse(row,col,val,sum(mb.^2),ncons);

c = zeros(sum(mb.^2), 1);
for i = 1:lsp
    ind = nbfind(sp, lsp, sp(:,i), n);
    c(loa{ind}) = coe(i)/length(loa{ind})*ones(length(loa{ind}),1);
end    
K.s = mb;
end