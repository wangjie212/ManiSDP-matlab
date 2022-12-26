function [cx, S, y, M] = spadmm(blk, At, C, b, dA, S, M)
gama = 1.99;
sigma = 1e-3;
maxiter = 5000;
mb = blk{1,2};
vmb = mb*(mb+1)/2;
A = [At{1}; At{2}];
dA = 1./dA;
iD = sparse(1:length(dA),1:length(dA),dA);
iA = iD - iD*At{2}'*(sparse(1:size(At{2},1),1:size(At{2},1),ones(size(At{2},1),1))+At{2}*iD* At{2}')^(-1)* At{2}*iD;
if isempty(S)
    S = zeros(size(A,1),1);
end
if isempty(M)
    M = zeros(size(A,1),1);
end
C = [mattovec(C{1}); full(C{2})];
i = 1;
error = 1;
while i <= maxiter && error >= 1e-6
    temp = 1/sigma*S + M - C;
    y = iA*(1/sigma*b - A'*temp);
    temp = S + sigma*(A*y - C);
    [V,D] = eig(vectomat(temp(1:vmb), mb), 'vector');
    psd = V*diag(max(0,D))*V';
%     [V,D,~] = poseig(vectomat(temp(1:vmb), mb));
%     psd = sparse(mb,mb);
%     for j = 1:length(find(D))
%         psd = psd + D(j)*V(:,j)*V(:,j)';
%     end
    ntemp = temp;
    ntemp(1:vmb) = mattovec(psd);
    M = 1/sigma*(ntemp - temp);
    temp = 1/sigma*S + M - C;
    y = iA*(1/sigma*b - A'*temp);
    S = S + gama*sigma*(M + A*y - C);
    if mod(i, 200) == 0
        cx = - C'*S;
        error = max([norm(A'*S-b)/(1+norm(b)), abs(C'*S-b'*y)/(1+abs(C'*S)+abs(b'*y)), norm(M+A*y-C)/(1+norm(C))]);
        disp(['spADMM iteration ' num2str(i) ': fval = ' num2str(cx,10) ', error = ' num2str(error,10)]);
    end
    i = i + 1;
end
end