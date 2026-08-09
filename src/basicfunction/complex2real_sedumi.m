%COMPLEX2REAL_SEDUMI Convert a complex SDP to a real SDP in SeDuMi format.
%
%   [AR,BR,CR,KR] = ...
%       COMPLEX2REAL_SEDUMI(A,B,C,K) converts
%
%       minimize    real(C'*X)
%       subject to  A*X = B,  X in K,
%
% from SeDuMi's complex/Hermitian representation to an equivalent problem
% containing real data and real cones only.  The conversion uses the
% economical realification in Theorem 2.2 of
%
%       J. Wang, "A More Efficient Reformulation of Complex SDP as Real
%       SDP", arXiv:2307.11599.
%
% For every n-by-n Hermitian PSD block H, the output has one unrestricted
% 2n-by-2n real symmetric PSD block
%
%              X = [X11, X12; X12', X22] >= 0,
%
% and the corresponding complex matrix is recovered by
%
%       H = (X11 + X22) + 1i*(X12' - X12).                 (1)
%
% Equivalently, if X21 = X12', then
%
%       H = (X11 + X22) + 1i*(X21 - X12).
%
% Unlike the standard structured realification, this formulation imposes
% no additional affine equations on X.
%
% INPUT
%   A,B,C,K are standard SeDuMi data.  In particular:
%     * K.scomplex lists the indices of the Hermitian PSD blocks;
%     * K.ycomplex lists rows for which A(row,:)*X = B(row) is a full
%       complex equality.  A row with imag(B(row)) ~= 0 is also treated as
%       complex, consistently with SeDuMi;
%     * all primal variables outside K.scomplex must be real.  Therefore
%       K.xcomplex must be absent or empty.
%   DAAT is diag(A*A'), as returned by cqpsos.m.
%
% OUTPUT
%   AR,BR,CR,KR are real SeDuMi data.  DAATR is diag(AR*AR').
%   Solve the converted problem with
%
%       [xreal,yreal,info] = sedumi(AR,BR,CR,KR);
%
%   MAP contains the sparse recovery map and indexing information:
%
%       xcomplex = MAP.T*xreal;
%
%   recovers the original primal vector, including every Hermitian block.
%   The original dual multiplier is recovered by
%
%       ycomplex = yreal(1:MAP.m_original);
%       q = MAP.complex_rows;
%       ycomplex(q) = ycomplex(q) + ...
%           1i*yreal(MAP.m_original + (1:numel(q)));
%
%   The first MAP.m_original equations in the real problem are the real
%   parts of all original equations; the remaining equations are the
%   imaginary parts of MAP.complex_rows, in that order.
%
% EXAMPLE: convert data produced by cqpsos.m
%   [A,b,c,K,dAAt] = cqpsos(Q,[],size(Q,1));
%   [Ar,br,cr,Kr,dAAt_r,map] = ...
%       complex2real_sedumi(A,b,c,K,dAAt);
%   [xr,yr,info] = sedumi(Ar,br,cr,Kr);
%   xc = map.T*xr;
%
% LIMITATION
%   The paper's construction concerns Hermitian PSD variables.  Complex
%   free, linear, Lorentz, or rotated-Lorentz variables (K.xcomplex) are
%   deliberately rejected rather than silently converted incorrectly.

function [Ar,br,cr,Kr] = complex2real_sedumi(A,b,c,K)

if nargin ~= 4
    error('complex2real_sedumi:InvalidInputCount', ...
        'Four inputs A, b, c, and K are required.');
end
if ~isstruct(K)
    error('complex2real_sedumi:InvalidCone', ...
        'K must be a SeDuMi cone structure.');
end
if ~isnumeric(A) || ndims(A) ~= 2
    error('complex2real_sedumi:InvalidA', ...
        'A must be a numeric matrix.');
end

[m,nA] = size(A);
b = b(:);
c = c(:);
if numel(b) ~= m
    error('complex2real_sedumi:DimensionMismatch', ...
        'numel(b) must equal size(A,1).');
end
% if ~isnumeric(dAAt) || ~isreal(dAAt) || numel(dAAt) ~= m || ...
%         any(~isfinite(dAAt(:))) || any(dAAt(:) < 0)
%     error('complex2real_sedumi:InvaliddAAt', ...
%         ['dAAt must be a finite, nonnegative real vector with ', ...
%          'size(A,1) entries.']);
% end
% dAAt_original = full(dAAt(:));

% Read the standard SeDuMi cone fields, treating absent fields as empty.
kf = cone_scalar(K,'f');
kl = cone_scalar(K,'l');
kq = cone_vector(K,'q');
kr = cone_vector(K,'r');
ks = cone_vector(K,'s');

n_prefix = kf + kl + sum(kq) + sum(kr);
n_original = n_prefix + sum(ks.^2);
if nA ~= n_original || numel(c) ~= n_original
    error('complex2real_sedumi:ConeDimensionMismatch', ...
        ['The column dimension implied by K is %d, whereas size(A,2) ', ...
         'is %d and numel(c) is %d.'],n_original,nA,numel(c));
end

if isfield(K,'xcomplex') && ~isempty(K.xcomplex)
    error('complex2real_sedumi:UnsupportedComplexVariables', ...
        ['Nonempty K.xcomplex is not supported.  This converter applies ', ...
         'the paper''s realification to Hermitian PSD blocks only.']);
end

ns = numel(ks);
if isfield(K,'scomplex')
    scomplex = normalize_indices(K.scomplex,ns,'K.scomplex');
else
    scomplex = zeros(0,1);
end
is_complex_psd = false(ns,1);
is_complex_psd(scomplex) = true;

if isfield(K,'ycomplex')
    ycomplex = normalize_indices(K.ycomplex,m,'K.ycomplex');
else
    ycomplex = zeros(0,1);
end
% This is SeDuMi's documented convention: a nonreal right-hand side also
% turns the corresponding row into a full complex equality.
complex_rows = unique([ycomplex; find(imag(b) ~= 0)],'sorted');

%% Construct x_original = T*x_real block by block.
% Variables preceding the PSD blocks, and already-real PSD blocks, pass
% through unchanged.  Each Hermitian n-by-n block receives the map (1).
Tblocks = cell(ns+1,1);
Tblocks{1} = speye(n_prefix);

old_psd_ranges = zeros(ns,2);
new_psd_ranges = zeros(ns,2);
old_cursor = n_prefix;
new_cursor = n_prefix;
ks_real = ks;

for j = 1:ns
    n = ks(j);
    old_psd_ranges(j,:) = old_cursor + [1,n^2];
    if is_complex_psd(j)
        Tblocks{j+1} = hermitian_recovery_map(n);
        ks_real(j) = 2*n;
    else
        Tblocks{j+1} = speye(n^2);
    end
    new_psd_ranges(j,:) = new_cursor + [1,ks_real(j)^2];
    old_cursor = old_cursor + n^2;
    new_cursor = new_cursor + ks_real(j)^2;
end

T = blkdiag(Tblocks{:});
n_real = size(T,2);

%% Transform the affine equations and objective.
% A uses no conjugation in A*x=b, whereas SeDuMi writes the objective as
% c'*x.  Consequently the two transformations are A*T and T'*c.
Awork = A*T;
cwork = T'*c;

Ar = [real(Awork); real_imag_part(Awork(complex_rows,:))];
br = [real(b); imag(b(complex_rows))];
cr = real(cwork);

% A SeDuMi PSD block represents a symmetric matrix.  Although the raw
% recovery map is equivalent on symmetric cone variables, Awork may encode
% a nonsymmetric coefficient matrix whose skew part vanishes identically.
% Remove that irrelevant skew part explicitly so that every PSD matrix
% encoded in every row of Ar (and in cr) is symmetric.
[Ar,cr] = symmetrize_psd_coefficients(Ar,cr,n_prefix,ks_real);

% Preserve sparsity whenever the input was sparse.
if issparse(A)
    Ar = sparse(Ar);
end
if issparse(c)
    cr = sparse(cr);
end

% Row squared norms used by diagonal normal-equation preconditioners.
% Since Ar is real, this is exactly diag(Ar*Ar') without forming Ar*Ar'.
% dAAt = full(sum(Ar.^2,2));

% Warn only if a genuinely contributing imaginary objective component is
% present.  Skew-symmetric imaginary coefficients on a symmetric PSD block
% vanish identically and are therefore harmless.
objective_imag_residual = objective_imaginary_residual(cwork,n_prefix,ks_real);
objective_scale = max(1,norm(cwork));
if objective_imag_residual > 1e-10*objective_scale
    warning('complex2real_sedumi:NonrealObjective', ...
        ['The objective is not real-valued on the stated real/Hermitian ', ...
         'cone to working precision; its imaginary part was discarded.']);
end

%% Output cone: every PSD block is now real symmetric.
Kr = K;
Kr.f = kf;
Kr.l = kl;
Kr.q = kq;
Kr.r = kr;
Kr.s = reshape(ks_real,size(ks));
remove = intersect(fieldnames(Kr),{'xcomplex','scomplex','ycomplex'});
if ~isempty(remove)
    Kr = rmfield(Kr,remove);
end

%% Recovery and audit information.
% map = struct();
% map.T = T;
% map.m_original = m;
% map.n_original = n_original;
% map.n_real = n_real;
% map.complex_rows = complex_rows;
% map.real_equation_rows = (1:m).';
% map.imag_equation_rows = m + (1:numel(complex_rows)).';
% map.s_original = ks;
% map.s_real = ks_real;
% map.scomplex_original = scomplex;
% map.old_psd_ranges = old_psd_ranges;
% map.new_psd_ranges = new_psd_ranges;
% map.dAAt_original = dAAt_original;
% map.dAAt_real = dAAt;

end

function T = hermitian_recovery_map(n)
% Return T such that vec(H)=T*vec(X), with X in S^(2n), and
% H=(X11+X22)+1i*(X21-X12).
N = 2*n;
[ii,jj] = ndgrid(1:n,1:n);
row = ii + (jj-1)*n;

col11 = ii     + (jj-1)*N;
col22 = (n+ii) + (n+jj-1)*N;
col21 = (n+ii) + (jj-1)*N;
col12 = ii     + (n+jj-1)*N;

T = sparse([row(:);row(:);row(:);row(:)], ...
           [col11(:);col22(:);col21(:);col12(:)], ...
           [ones(n^2,1);ones(n^2,1); ...
            1i*ones(n^2,1);-1i*ones(n^2,1)], ...
           n^2,N^2);
end

function Z = real_imag_part(Zcomplex)
% Keep this helper explicit so that an empty complex-row set produces a
% correctly sized real sparse matrix in both MATLAB and Octave.
Z = imag(Zcomplex);
end

function [A,c] = symmetrize_psd_coefficients(A,c,n_prefix,ks)
% Canonicalize all real SDP coefficient blocks.  If a row segment is
% vec(M)', replace it by vec((M+M')/2)'.  This leaves the linear functional
% unchanged on a symmetric PSD variable and makes the encoded matrix
% explicitly symmetric.
cursor = n_prefix;
for j = 1:numel(ks)
    n = ks(j);
    ind = cursor + (1:n^2);

    % With column-major vectorization, transpose_order maps vec(M) to
    % vec(M').
    transpose_order = reshape(1:n^2,n,n).';
    transpose_order = transpose_order(:);

    Ablock = A(:,ind);
    A(:,ind) = 0.5*(Ablock + Ablock(:,transpose_order));

    cblock = c(ind);
    c(ind) = 0.5*(cblock + cblock(transpose_order));

    cursor = cursor + n^2;
end
end

function value = cone_scalar(K,name)
if ~isfield(K,name) || isempty(K.(name))
    value = 0;
else
    value = K.(name);
end
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        value < 0 || value ~= floor(value)
    error('complex2real_sedumi:InvalidConeField', ...
        'K.%s must be a nonnegative integer scalar.',name);
end
end

function value = cone_vector(K,name)
if ~isfield(K,name) || isempty(K.(name))
    value = zeros(0,1);
else
    value = K.(name);
end
if ~isnumeric(value) || ~isreal(value) || any(value(:) <= 0) || ...
        any(value(:) ~= floor(value(:)))
    error('complex2real_sedumi:InvalidConeField', ...
        'Every entry of K.%s must be a positive integer.',name);
end
end

function ind = normalize_indices(raw,nmax,name)
if isempty(raw)
    ind = zeros(0,1);
    return
end
if islogical(raw)
    if numel(raw) ~= nmax
        error('complex2real_sedumi:InvalidConeField', ...
            'Logical %s must have %d entries.',name,nmax);
    end
    ind = find(raw(:));
else
    ind = raw(:);
    if ~isnumeric(ind) || ~isreal(ind) || any(ind ~= floor(ind)) || ...
            any(ind < 1) || any(ind > nmax)
        error('complex2real_sedumi:InvalidConeField', ...
            '%s contains an invalid index.',name);
    end
    ind = unique(ind,'sorted');
end
end

function residual = objective_imaginary_residual(c,n_prefix,ks)
% Measure only imaginary coefficients that can act on the real cone.
residual_squared = norm(imag(c(1:n_prefix)))^2;
cursor = n_prefix;
for j = 1:numel(ks)
    n = ks(j);
    ind = cursor + (1:n^2);
    M = reshape(imag(c(ind)),n,n);
    Msym = (M+M.')/2;
    residual_squared = residual_squared + norm(Msym,'fro')^2;
    cursor = cursor + n^2;
end
residual = sqrt(residual_squared);
end
