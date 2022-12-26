function [SDP,info] = dense_sdp_relax_binary(problem,kappa)

if ~isfield(problem,'vars'); error('Please provide variables of the POP.'); end
if ~isfield(problem,'objective'); error('Please provide objective function.'); end
if ~isfield(problem,'equality'); problem.equality = []; end

fprintf('\n======================= Dense SDP Relaxation =======================\n')

%% copy POP data and decide relaxation order
x       = problem.vars;
f       = problem.objective;
h       = problem.equality;

if kappa == 1; v = [1;x]; end
if kappa > 1; v = [1;x;monomials(x(end:-1:1),2:kappa)]; end

ind = true(length(v), 1);
for i = 1:length(v)
    mon = v(i);
    if deg(mon) == 2 && length(mon.var) == 1
        ind(i) = false;
    end
end
v = v(ind);

J = diff(f,x);
info.f                   = f;
info.J                   = J'; % for evaluate cost gradient
info.v                   = v; % for lifting
info.d                   = length(x);

xv = x.var;
%% compute multipliers
fprintf('Compute equality and inequality multipliers ... ')
kappa2  = 2*kappa;
pop     = [mykron(v,v);f];
% equalities
l_h     = length(h);
dim_h   = 0;
if l_h > 0
    for i = 1:l_h
        hi      = h(i);
        deghi   = deg(hi);
        lamhi   = monomials(x,0:(kappa2 - deghi));
        ind = true(length(lamhi), 1);
        for j = 1:length(lamhi)
            mon = lamhi(j);
            if sum(mon.var == xv(i)) == 1
                ind(j) = false;
            end
        end
        lamhi = lamhi(ind);
        pop     = [pop;hi*lamhi];
        dim_h   = dim_h + length(lamhi);
    end
end
fprintf('Done.\n')

%% decompose polynomials to generate SDP data
[~,degmat,coef_all] = decomp(pop);
coef_all            = coef_all';

n                   = length(v); % size of the moment matrix
ndelta              = triangle_number(n);
nterms              = size(degmat,1);
m_mom               = ndelta - nterms;
m_h                 = dim_h; % number of constraints due to equality
m                   = m_mom + m_h + 1;

coef_mom    = coef_all(:,1:n^2);
coef_mom    = coef_mom';

B           = {};
B_normalize = {};
A           = {};

fprintf('Build moment constraints ... ')
for i = 1:nterms
    if rem(i,10000) == 1
        fprintf('%d/%d ',i,nterms);
    end

    [row,~,~]   = find(coef_mom(:,i));
    SDP_coli    = floor((row-1)./n) + 1;
    SDP_rowi    = mod(row-1,n) + 1;
    nnz         = length(SDP_rowi);
    
    Bi          = sparse(SDP_rowi,SDP_coli,ones(nnz,1),n,n);
    B{end+1}    = Bi;
    B_normalize{end+1} = Bi/nnz;
    
    mask_triu   = (SDP_rowi >= SDP_coli);
    si          = SDP_rowi(mask_triu);
    sj          = SDP_coli(mask_triu);

    nnz_triu    = length(si);
    
    if nnz_triu > 1
        [~,base_idx]        = max(sj);
        si_base             = si(base_idx);
        sj_base             = sj(base_idx);
        
        si_nonbase          = si; 
        si_nonbase(base_idx)= [];
        sj_nonbase          = sj; 
        sj_nonbase(base_idx)= [];
        
        is_base_diag        = (si_base == sj_base);
        
        if is_base_diag
            A_si            = [si_base];
            A_sj            = [sj_base];
            A_v             = [1];
        else
            A_si            = [si_base,sj_base];
            A_sj            = [sj_base,si_base];
            A_v             = [0.5,0.5];
        end
        
        for nonbase_idx = 1:length(si_nonbase)
            is_nonbase_diag = (si_nonbase(nonbase_idx) == sj_nonbase(nonbase_idx));
            if is_nonbase_diag
                A_sii       = [A_si,si_nonbase(nonbase_idx)];
                A_sjj       = [A_sj,sj_nonbase(nonbase_idx)];
                A_vv        = [A_v,-1];
            else
                A_sii       = [A_si,si_nonbase(nonbase_idx),sj_nonbase(nonbase_idx)];
                A_sjj       = [A_sj,sj_nonbase(nonbase_idx),si_nonbase(nonbase_idx)];
                A_vv        = [A_v,-0.5,-0.5];
            end
            A_temp          = sparse(A_sii,A_sjj,A_vv,n,n);
            
            A{end+1}        = A_temp;
        end
    end
end
fprintf('Done.\n')

%% Build A's associated with localizing equality constraints
if m_h == 0
    % Do nothing
    A_local = {};
else
    coef_loc    = coef_all(:,1+n^2+(1:m_h));
    A_local     = {};
    fprintf('Build localizing equality constraints ... ')
    for i = 1:m_h
        if rem(i,10000) == 1
            fprintf('%d/%d ',i,m_h);
        end
        
        [rowi,~,vi] = find(coef_loc(:,i));
        
        Ai      = sparse(n,n);
        for j   = 1:length(rowi)
            Ai  = Ai + vi(j) * B_normalize{rowi(j)};
        end
        A_local{end+1} = Ai;
    end
end
fprintf('Done.\n')

%% Leading A
A0 = sparse([1],[1],[1],n,n);
% Combine all A for the main block
A = [{A0},A_local,A];

%% Now build the cost matrix
coef_cost   = coef_all(:,n^2+1);
[row,~,v]   = find(coef_cost);
C           = sparse(n,n);
fprintf('Build cost matrix C ... ')
for i = 1:length(row)
    if rem(i,1000) == 1
        fprintf('%d/%d ',i,length(row));
    end
    C       = C + v(i) * B_normalize{row(i)};
end
fprintf('Done.\n')

%% convert to SDPT3 format
fprintf('Generate SDPT3 data ... ')
blk{1,1}        = 's';
blk{1,2}        = n;
b               = sparse([1],[1],[1],m,1);
At0             = sparsesvec(blk(1,:),A);
At              = {At0};
C       = {C};

SDP.blk = blk;
SDP.At  = At;
SDP.C   = C;
SDP.b   = b;

fprintf('Done.\n')

%% convert to sedumi format
fprintf('Generate sedumi data ... ')
sK.s  = [n];

A0t     = sparsevec(blk(1,:),A);
At      = {A0t};

sdata.K     = sK;
sdata.At    = cat(1,At{:});
sdata.b     = b;

sc          = [];
for i = 1:length(SDP.C)
    sc      = [sc;sparsevec(blk(i,:),SDP.C(i))];
end
sdata.c     = sc;

SDP.sedumi   = sdata;

fprintf('Done.\n')
fprintf('====================================================================\n\n\n')
info.kappa  = kappa;

end