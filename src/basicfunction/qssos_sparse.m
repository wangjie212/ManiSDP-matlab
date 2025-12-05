% Generate the second-order SOS-SDP relaxation for the sparse quartic sphere program:
% inf  c1'*[x_{I_1}]_4 + ... + ct'*[x_{I_t}]_4
% s.t. x_{I_1}^2 = 1, ..., x_{I_t}^2 = 1.
% Output sedumi format data.
% dAAt: the diagnal of A*A'

function [A, b, c, K, dAAt] = qssos_sparse(n, I, coe)
t = length(I);
mb = zeros(t, 1);
sp = [];
cql = zeros(t, 1);
for i = 1:t
    basis{i} = get_basis(n, 2, I{i});
    mb(i) = size(basis{i}, 2);
    sp = [sp get_basis(n, 4, I{i})];
    cql(i) = length(I{i});
end
sp = unique(sp', 'rows');
sp = sortrows(sp)';
lsp = size(sp, 2);
nz = sum(mb.^2 + mb.*(cql+1)) + 1;
row = zeros(nz, 1);
col = zeros(nz, 1);
val = ones(nz, 1);
dAAt = zeros(lsp, 1);
row(1) = 1;
col(1) = 1;
ind = 2;
for k = 1:t
    ss = sum(mb(1:k-1)) + 1;
    for i = 1:mb(k)
        for j = 1:cql(k)+1
            temp = zeros(n, 1);
            if j < cql(k)+1
               temp(I{k}(j)) = 2;
            else
               val(ind) = -1;
            end
            bi = basis{k}(:,i) + temp;
            locb = nbfind(sp, lsp, bi, n);
            row(ind) = locb;
            col(ind) = ss + i;
            ind = ind + 1;
        end
    end
end
for k = 1:t
    ss = sum(mb(1:k-1).^2) + sum(mb) + 1;
    for i = 1:mb(k)
        for j = i:mb(k)
            bi = basis{k}(:,i) + basis{k}(:,j);
            locb = nbfind(sp, lsp, bi, n);
            row(ind) = locb;
            col(ind) = ss + (i-1)*mb(k) + j;
            if j > i
                row(ind+1) = locb;
                col(ind+1) = ss + (j-1)*mb(k) + i;
                dAAt(locb) = dAAt(locb) + 2;
                ind  = ind + 2;
            else
                dAAt(locb) = dAAt(locb) + 1;
                ind = ind + 1;
            end
        end
    end
end
nvar = sum(mb.^2 + mb) + 1;
A = sparse(row,col,val,lsp,nvar);
b = coe;
K.f = sum(mb) + 1;
K.s = mb';
c = zeros(nvar, 1);
c(1) = 1;
end
