% Generate the second-order sparse moment-SDP relaxation for the sensor network localization problem.
% Output sedumi format data.

function [At, b, c, K] = snl_mom_sparse(f, x, cliques)
n = length(x);
[~, tsupp, tcoe] = decomp(f);
[~, loc, coe] = find(tcoe);
supp = tsupp(loc, :)';
t = length(cliques);
mb = zeros(t, 1);
sp = [];
for i = 1:t
    basis{i} = get_basis(n, 2, cliques{i});
    mb(i) = size(basis{i}, 2);
    sp = [sp get_basis(n, 4, cliques{i})];
end
sp = unique(sp', 'rows');
sp = sortrows(sp)';
lsp = size(sp, 2);
mm = cell(lsp,1);
for k = 1:t
    for i = 1:mb(k)
        for j = i:mb(k)
            bi = basis{k}(:,i) + basis{k}(:,j);
            ind = nbfind(sp, lsp, bi, n);
            mm{ind} = [mm{ind} [i;j;k]];
        end
    end
end
ncons = sum(mb.*(mb+1)/2) - lsp + sum(mb) + 1;
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
          temp = sum(mb(1:mm{i}(3,j)-1).^2);
          loa{i}(2*j-1:2*j) = [temp+(mm{i}(2,j)-1)*mb(mm{i}(3,j)) + mm{i}(1,j);temp+(mm{i}(1,j)-1)*mb(mm{i}(3,j)) + mm{i}(2,j)];
     end
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
At = sparse(row,col,val,sum(mb.^2),ncons);

c = zeros(sum(mb.^2), 1);
for i = 1:size(supp,2)
    ss = [];
    ind = nbfind(sp, lsp, supp(:,i), n);
    for k = 1:size(mm{ind},2)
        if mm{ind}(1,k) == mm{ind}(2,k)
            ss = [ss; loa{ind}(2*k)];
        else
            ss = [ss; loa{ind}(2*k-1:2*k)];
        end
    end
    c(ss) = coe(i)/length(ss)*ones(length(ss),1);
end
K.s = mb';
end
