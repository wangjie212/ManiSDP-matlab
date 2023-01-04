% Generate the second-order moment-SDP relaxation for the binary quadratic program:
% Min  x'*Q*x + x'*e
% s.t. x_i^2 = 1, i = 1,...,n.
% Output sedumi format data.
% If unitdiag = 1, then include unit diagonal constraints.

function [At, b, c, mb] = bqpmom(n, Q, e, unitdiag)
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
ncons = mb*(mb+1)/2 - lsp - mb + n*(mb-n-1);
if unitdiag == 1
    ncons = ncons + mb;
end
row = [];
col = [];
val = [];
b = sparse(ncons, 1);
l = 1;
if unitdiag == 1
    b(1:mb) = ones(mb,1);
    for i = 1:mb
        row = [row; (i-1)*mb+i];
        col = [col; i];
        val = [val; 1];
    end
    l = l + mb;
end
loa = cell(lsp,1);
lob = cell(lsp,1);
for i = 1:lsp
     loca0 = (mm{i}(2,1)-1)*mb + mm{i}(1,1);
     locb0 = (mm{i}(1,1)-1)*mb + mm{i}(2,1);
     loa{i} = [loca0];
     lob{i} = [locb0];
     for j = 2:size(mm{i},2)
          loca = (mm{i}(2,j)-1)*mb + mm{i}(1,j);
          locb = (mm{i}(1,j)-1)*mb + mm{i}(2,j);
          loa{i} = [loa{i};loca];
          lob{i} = [lob{i};locb];
          row = [row; loca0; locb0; loca; locb];
          col = [col; l; l; l; l];
          val = [val; 0.5; 0.5; -0.5; -0.5];
          l = l + 1;
     end
end
for k = 1:n
    for i = 2:mb
        if basis(k, i) == 0
            bi = basis(:,i);
            bi(k) = 2;
            ind1 = bfind(sp, lsp, bi, n);
            ind2 = bfind(sp, lsp, basis(:,i), n);
            row = [row; loa{ind1}; lob{ind1}; loa{ind2}; lob{ind2}];
            col = [col; l*ones(2*(length(loa{ind1})+length(loa{ind2})),1)];
            val = [val; 1/(2*length(loa{ind1}))*ones(2*length(loa{ind1}),1); -1/(2*length(loa{ind2}))*ones(2*length(loa{ind2}),1)];
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