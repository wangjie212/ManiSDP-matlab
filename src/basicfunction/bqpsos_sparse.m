% Generate the second-order SOS-SDP relaxation for the binary sparse quadratic program:
% inf  c1'*[x_{I_1}]_2 + ... + ct'*[x_{I_t}]_2
% s.t. x_i^2 = 1, i = 1,...,n.
% Output sedumi format data.
% dAAt: the diagnal of A*A'

function [A, b, c, K, dAAt] = bqpsos_sparse(n, I, coe)
t = length(I);
mb = zeros(t, 1);
sp = [];
for i = 1:t
    basis{i} = get_basis(n, 2, I{i});
    ind = true(size(basis{i}, 2), 1);
    ind(sum(basis{i}>1)> 0) = false;
    basis{i} = basis{i}(:,ind);
    mb(i) = size(basis{i}, 2);
    temp = get_basis(n, 4, I{i});
    ind = true(size(temp, 2), 1);
    ind(sum(temp>1)> 0) = false;
    sp = [sp temp(:,ind)];
end
sp = unique(sp', 'rows');
sp = sortrows(sp)';
lsp = size(sp, 2);
nz = sum(mb.^2) + 1;
row = ones(nz, 1);
col = zeros(nz, 1);
val = ones(nz, 1);
dAAt = zeros(lsp, 1);
dAAt(1) = sum(mb);
col(1) = 1;
ind = 2;
for k = 1:t
    ss = sum(mb(1:k-1).^2) + 1;
    col(ind:ind+mb(k)-1) = ss + (0:(mb(k)-1))*mb(k)+(1:mb(k));
    ind = ind + mb(k);
    for i = 1:mb(k)
        for j = i+1:mb(k)
            bi = mod(basis{k}(:,i) + basis{k}(:,j), 2);
            locb = nbfind(sp, lsp, bi, n);
            row(ind) = locb;
            row(ind+1) = locb;
            col(ind) = ss + (i-1)*mb(k)+j;
            col(ind+1) = ss + (j-1)*mb(k)+i;
            dAAt(locb) = dAAt(locb) + 2;
            ind = ind + 2;
        end
    end
end
A = sparse(row,col,val,lsp,nz);
b = coe;
K.f = 1;
K.s = mb';
c = zeros(nz, 1);
c(1) = 1;
end
