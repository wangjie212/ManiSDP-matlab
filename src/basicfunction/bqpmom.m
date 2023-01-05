% Generate the second-order moment-SDP relaxation for the binary quadratic program:
% Min  x'*Q*x + x'*e
% s.t. x_i^2 = 1, i = 1,...,n.
% Output sedumi format data.

function [At, b, c, mb] = bqpmom(n, Q, e)
basis = get_basis(n, 2);
ind = true(size(basis, 2), 1);
for i = 1:size(basis, 2)
    if sum(basis(:,i) > 1) > 0
        ind(i) = false;
    end
end
basis = basis(:,ind);
mb = size(basis, 2);
sp = get_basis(n, 4);
ind = true(size(sp, 2), 1);
for i = 1:size(sp, 2)
    if sum(sp(:,i) > 2) > 0 || sum(mod(sp(:,i),2)) == 0
        ind(i) = false;
    end
end
sp = sp(:,ind);
lsp = size(sp, 2);
mm = cell(lsp,1);
for i = 1:mb
    for j = i+1:mb
         bi = basis(:,i) + basis(:,j);
         ind = bfind(sp, lsp, bi, n);
         mm{ind} = [mm{ind} [i;j]];
    end
end
ncons = mb*(mb+1)/2 - lsp + n*(mb-1) - mb + 1;
row = [1];
col = [1];
val = [1];
b = sparse(ncons, 1);
b(1) = 1;
for i = 2:n+1
    row = [row; 1; (i-1)*mb+i];
    col = [col; i; i];
    val = [val; 1; -1];
end
l = n + 2;
ind = [1:n];
for i = n+2:mb
    cc = ind(basis(:,i)==1) + 1;
    row = [row; (cc(1)-1)*mb+cc(1); (i-1)*mb+i; (cc(2)-1)*mb+cc(2); (i-1)*mb+i];
    col = [col; l; l; l+1; l+1];
    val = [val; 1; -1; 1; -1];
    l = l + 2;
end
loa = cell(lsp,1);
for i = 1:lsp
     loa{i} = zeros(2*size(mm{i},2),1);
     for j = 1:size(mm{i},2)
          loa{i}(2*j-1:2*j) = [(mm{i}(2,j)-1)*mb + mm{i}(1,j);(mm{i}(1,j)-1)*mb + mm{i}(2,j)];
     end
end
for k = 1:n
    for i = 2:size(basis, 2)
        if basis(k, i) == 0
            bi = basis(:,i);
            bi(k) = 2;
            ind1 = bfind(sp, lsp, bi, n);
            ind2 = bfind(sp, lsp, basis(:,i), n);
            row = [row; loa{ind1}; loa{ind2}];
            col = [col; l*ones(length(loa{ind1})+length(loa{ind2}),1)];
            val = [val; 1/length(loa{ind1})*ones(length(loa{ind1}),1); -1/length(loa{ind2})*ones(length(loa{ind2}),1)];
            l = l + 1;          
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
At = sparse(row,col,val,mb^2,ncons);

row = [2:n+1]';
col = [2:n+1]';
val = diag(Q);
for i = 1:n
    for j = 1:size(mm{i},2)
        row = [row; mm{i}(1,j); mm{i}(2,j)];
        col = [col; mm{i}(2,j); mm{i}(1,j)];
    end
    val = [val; e(i)/(2*size(mm{i},2))*ones(2*size(mm{i},2),1)];
end
ind = n + 1;
for i = 2:n
    for j = 1:i-1
        for k = 1:size(mm{ind},2)
            row = [row; mm{ind}(1,k); mm{ind}(2,k)];
            col = [col; mm{ind}(2,k); mm{ind}(1,k)];
        end
        val = [val; Q(j,i)/size(mm{ind},2)*ones(2*size(mm{ind},2),1)];
        ind = ind + 1;
    end
end    
C = sparse(row,col,val,mb,mb);
c = C(:);
end