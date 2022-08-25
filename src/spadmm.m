function [X, y, S] = spadmm(blk, At, C, b, dA)
gama = 1.99;
sigma = 0.001;
maxiter = 5000;
mb = blk{1,2};
vmb = mb*(mb+1)/2;
A = [At{1}; At{2}];
dA = 1./dA;
iD = sparse(1:length(dA),1:length(dA),dA);
iA = iD - iD*At{2}'*(sparse(1:size(At{2},1),1:size(At{2},1),ones(size(At{2},1),1))+At{2}*iD* At{2}')^(-1)* At{2}*iD;
X = sparse(size(A,1),1);
S = sparse(size(A,1),1);
C = sparse([mattovec(C{1}); C{2}]);
i = 1;
error = 1;
while i <= maxiter && error >= 0.1
    temp = 1/sigma*X + S - C;
    y = iA*(1/sigma*b - A'*temp);
    temp = X + sigma*(A*y - C);
    [V,D] = eig(vectomat(temp(1:vmb), mb), 'vector');
    psd = V*diag(max(0,D))*V';
%     [V,D,~] = poseig(vectomat(temp(1:vmb), mb));
%     psd = sparse(mb,mb);
%     for j = 1:length(find(D))
%         psd = psd + D(j)*V(:,j)*V(:,j)';
%     end
    ntemp = temp;
    ntemp(1:vmb) = mattovec(psd);
    S = 1/sigma*(ntemp - temp);
    temp = 1/sigma*X + S - C;
    y = iA*(1/sigma*b - A'*temp);
    X = X + gama*sigma*(S + A*y - C);
    if mod(i, 100) == 0
        cx = - C'*X;
        error = max(norm(X'*S), norm(S+A*y-C));
        disp(['spADMM iteration ' num2str(i) ': fval = ' num2str(cx,10) ', error = ' num2str(error,10)]);
    end
    i = i + 1;
end
end