function [At,b,c,K] = cquarticmom_nounitdiag(n,Q)
%CQUARTICMOM_REDUCED_NOUNITDIAG Reduced order-2 complex moment SDP.
%
%   [At,b,c,K] = cquarticmom_reduced_nounitdiag(n,Q)
%
% constructs the affine moment-consistency equations for
%
%       minimize    v_2(z)'*Q*v_2(z)
%       subject to  abs(z_i)^2 = 1,  i=1,...,n,
%
% using the reduced full mixed basis
%
%       {z^alpha*conj(z)^beta : |alpha|+|beta| <= 2}
%           \ {z_i*conj(z_i) : i=1,...,n}.
%
% Here v_2(z)={z^gamma: |gamma|<=2}.  The Hermitian moment block has
%
%       mb = nchoosek(2*n+2,2)-n
%
% rows and columns.
%
% IMPORTANT
% ---------
% The equations diag(X)=1 are deliberately NOT included in At and b.
% They must be enforced separately, for example by ManiDSDP_unitdiag.
% Consequently, the data returned by this function do not by themselves
% describe the original relaxation when passed directly to SeDuMi.
%
% To restore a standard SeDuMi model, append the unit-diagonal equations:
%
%   didx   = sub2ind([mb,mb],1:mb,1:mb);
%   AtDiag = sparse(didx,1:mb,1,mb^2,mb);
%   AtSed  = [At,AtDiag];
%   bSed   = [b;ones(mb,1)];
%   [x,y,info] = sedumi(AtSed,bSed,c,K);
%
% SeDuMi convention:
%
%       minimize    c'*x
%       subject to  At'*x=b,  x in K,
%
% where x=vec(X) and X is Hermitian positive semidefinite.
%
% Inputs
%   n       positive integer.
%   Q       mh-by-mh Hermitian matrix, mh=nchoosek(n+2,2), indexed by
%           1,z_i,z_i*z_j (1<=i<=j<=n).
%
% Outputs
%   At      affine moment-consistency operator; no unit-diagonal columns.
%   b       zero right-hand side.
%   c       vec(C), the objective vector.
%   K       one mb-by-mb complex Hermitian PSD SeDuMi block.
%   mb      reduced mixed-basis size.
%   alpha   holomorphic exponents of the reduced mixed basis.
%   beta    antiholomorphic exponents of the reduced mixed basis.
%   holind  positions of the holomorphic basis v_2(z).
%   dsp     retained Laurent-exponent representatives.
%   C       Hermitian objective matrix.
%
% Example with ManiDSDP_unitdiag (adapt the call to the local interface):
%
%   n  = 3;
%   mh = nchoosek(n+2,2);
%   R  = randn(mh)+1i*randn(mh);
%   Q  = (R+R')/2;
%   [At,b,c,K] = cquarticmom_nounitdiag(n,Q);
%   A = At';                       % row-wise affine operator if required
%   % diag(X)=1 is handled structurally by ManiDSDP_unitdiag.

if nargin < 2
    error('cquarticmom_reduced_nounitdiag:NotEnoughInputs', ...
        'Both n and Q are required.');
end
if ~isscalar(n) || ~isnumeric(n) || ~isreal(n) || ...
        n < 1 || n ~= floor(n)
    error('cquarticmom_reduced_nounitdiag:InvalidDimension', ...
        'n must be a positive integer.');
end
if ~isnumeric(Q) || ndims(Q) ~= 2
    error('cquarticmom_reduced_nounitdiag:InvalidQ', ...
        'Q must be a numeric matrix.');
end

mh = nchoosek(n+2,2);
if ~isequal(size(Q),[mh,mh])
    error('cquarticmom_reduced_nounitdiag:InvalidQSize', ...
        'Q must be %d-by-%d because length(v_2(z))=%d.',mh,mh,mh);
end

hermerr = norm(Q-Q','fro');
qscale = max(1,norm(Q,'fro'));
if hermerr > 1e-12*qscale
    error('cquarticmom_reduced_nounitdiag:NonHermitianQ', ...
        'Q must be Hermitian.');
end
Q = (Q+Q')/2;

%% Reduced full mixed basis of total degree at most two
[alpha,beta] = mixed_basis_order_two_nounitdiag(n);

% Remove precisely z_i*conj(z_i), i=1,...,n.
is_modulus_square = all(alpha==beta,1) & (sum(alpha,1)==1);
if nnz(is_modulus_square) ~= n
    error('cquarticmom_reduced_nounitdiag:InternalBasisError', ...
        'Could not identify all z_i*conj(z_i) monomials.');
end
alpha(:,is_modulus_square) = [];
beta(:,is_modulus_square) = [];

mb = size(alpha,2);
if mb ~= nchoosek(2*n+2,2)-n
    error('cquarticmom_reduced_nounitdiag:InternalBasisError', ...
        'Unexpected reduced mixed-basis size.');
end

degs = alpha-beta;
holind = find(all(beta==0,1));
if numel(holind) ~= mh
    error('cquarticmom_reduced_nounitdiag:InternalBasisError', ...
        'Could not identify the holomorphic degree-two subbasis.');
end

%% Group moment-matrix entries by Laurent exponent difference
nentry = mb^2;
entry_delta = zeros(nentry,n);
for k = 1:n
    delta_k = degs(k,:).'-degs(k,:);
    entry_delta(:,k) = delta_k(:);
end

[all_delta,~,entry_group] = unique(entry_delta,'rows');
ngroup = size(all_delta,1);
group_count = accumarray(entry_group,1,[ngroup,1]);

[has_opposite,opposite_group] = ismember(-all_delta,all_delta,'rows');
if ~all(has_opposite) || any(opposite_group(opposite_group)~=(1:ngroup)')
    error('cquarticmom_reduced_nounitdiag:InternalGroupError', ...
        'Laurent-exponent groups are not closed under conjugation.');
end

%% Objective matrix
% Embed Q in the holomorphic block and average over entries representing
% the same Laurent moment.
C0 = complex(zeros(mb,mb));
C0(holind,holind) = Q;
C0vec = C0(:);

group_sum = complex(zeros(ngroup,1));
for ell = 1:nentry
    g = entry_group(ell);
    group_sum(g) = group_sum(g)+C0vec(ell);
end
C = reshape(group_sum(entry_group)./group_count(entry_group),mb,mb);
C = (C+C')/2;
c = C(:);

%% Moment-consistency equations, excluding diag(X)=1
maxcons = nentry;
triplet_row = zeros(2*maxcons,1);
triplet_col = zeros(2*maxcons,1);
triplet_val = zeros(2*maxcons,1);
is_complex_equality = false(maxcons,1);

ncons = 0;
nnzAt = 0;

[~,group_order] = sort(entry_group);
group_start = [1;find(diff(entry_group(group_order)))+1];
group_end = [group_start(2:end)-1;nentry];

dsp = zeros(1,n);

for g = 1:ngroup
    delta = all_delta(g,:);
    first = find(delta~=0,1,'first');
    members = group_order(group_start(g):group_end(g));

    if isempty(first)
        % Keep only the upper triangle of the zero-Laurent group.
        % All diagonal members are skipped: diag(X)=1 is imposed by the
        % structural unit-diagonal solver.  Off-diagonal zero moments are
        % equated to X(1,1), which the structural constraint fixes to one.
        aind = mod(members-1,mb)+1;
        qind = floor((members-1)/mb)+1;
        members = members(aind<=qind);

        rep = 1;  % X(1,1), fixed externally to one
        if ~any(members==rep)
            error('cquarticmom_reduced_nounitdiag:InternalGroupError', ...
                'The zero group does not contain X(1,1).');
        end

        for t = 1:numel(members)
            current = members(t);
            ai = mod(current-1,mb)+1;
            qi = floor((current-1)/mb)+1;

            if ai==qi
                continue;             % remove every unit-diagonal equation
            end

            [ncons,nnzAt,triplet_row,triplet_col,triplet_val] = ...
                append_equality_nounitdiag(ncons,nnzAt,triplet_row, ...
                triplet_col,triplet_val,current,rep);
            is_complex_equality(ncons) = true;
        end

    elseif delta(first)>0
        % Retain one of each conjugate pair {delta,-delta}.
        dsp(end+1,:) = delta; %#ok<AGROW>
        rep = members(1);

        for t = 2:numel(members)
            current = members(t);
            [ncons,nnzAt,triplet_row,triplet_col,triplet_val] = ...
                append_equality_nounitdiag(ncons,nnzAt,triplet_row, ...
                triplet_col,triplet_val,current,rep);
            is_complex_equality(ncons) = true;
        end
    end
end

At = sparse(triplet_row(1:nnzAt),triplet_col(1:nnzAt), ...
    triplet_val(1:nnzAt),nentry,ncons);
b = zeros(ncons,1);

%% SeDuMi cone description
K = struct();
K.f = 0;
K.l = 0;
K.q = [];
K.r = [];
K.s = mb;
K.scomplex = 1;
K.ycomplex = find(is_complex_equality(1:ncons));

end

function [alpha,beta] = mixed_basis_order_two_nounitdiag(n)
% Exponents of 1,w_i,w_i*w_j, where w=[z;conj(z)].
mbfull = nchoosek(2*n+2,2);
alpha = zeros(n,mbfull);
beta = zeros(n,mbfull);
p = 1;

for i = 1:2*n
    p = p+1;
    if i<=n
        alpha(i,p) = 1;
    else
        beta(i-n,p) = 1;
    end
end

for i = 1:2*n
    for j = i:2*n
        p = p+1;
        if i<=n
            alpha(i,p) = alpha(i,p)+1;
        else
            beta(i-n,p) = beta(i-n,p)+1;
        end
        if j<=n
            alpha(j,p) = alpha(j,p)+1;
        else
            beta(j-n,p) = beta(j-n,p)+1;
        end
    end
end

if p~=mbfull || any(sum(alpha+beta,1)>2)
    error('cquarticmom_reduced_nounitdiag:InternalBasisError', ...
        'Unexpected full mixed basis.');
end
end

function [ncons,nnzAt,row,col,val] = append_equality_nounitdiag( ...
    ncons,nnzAt,row,col,val,current,representative)
% Add X(current)-X(representative)=0 to At'*vec(X)=b.
ncons = ncons+1;
row(nnzAt+1:nnzAt+2) = [current;representative];
col(nnzAt+1:nnzAt+2) = ncons;
val(nnzAt+1:nnzAt+2) = [1;-1];
nnzAt = nnzAt+2;
end
