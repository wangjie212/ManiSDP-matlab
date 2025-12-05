# ManiSDP
ManiSDP (v0.3.0) aims to solve the following low-rank semidefinite program (SDP) via manifold optimization:

$$\inf_{X\succeq0}\langle C, X\rangle{\ }{\ }\text{s.t.}{\ }{\ }\mathcal{A}(X)=b,{\ }\mathcal{B}(X)=d,$$

where the linear constraints $\mathcal{A}(X)=b$ are arbitrary while the linear constraints $\mathcal{B}(X)=d$, if present, are assumed to define certain manifold structure. Here, low-rank means the SDP admits a low-rank optimal solution.

## Dependencies
- [MATLAB](https://ww2.mathworks.cn/products/matlab.html?s_tid=hp_products_matlab)
- [Manopt](https://github.com/NicolasBoumal/manopt)

## Usage
**Note:** The optimal setting of parameters in ManiSDP is highly problem-dependent. We encourage the users to find the optimal parameters (defined in `options`) via preliminary experiments before running large-scale cases. In our experiences, the parameters (`sigma0`, `alpha`, `theta`, `TR_maxiter`, `TR_maxinner`, `tau1`, `tau2`) have a significant influence on the performance of ManiSDP.  

ManiSDP accepts [SeDuMi](https://sedumi.ie.lehigh.edu/) format data.  

So far, ManiSDP supports four types of SDPs.

### SDPs with only unit diagonal constraints
$$\inf_{X\succeq0}\langle C, X\rangle{\ }{\ }\text{s.t.}{\ }{\ }X_{ii}=1,{\ }\text{for}{\ }i=1,\ldots,n.$$

```matlab
clear options;
options.tol = 1e-8;
[sol, opt, data] = ManiSDP_onlyunitdiag(C, options);
```

**Input**  
`C`: cost matrix  
`options`:  
\- `tol` (=1e-8 by default): tolerance of maximal KKT residues  
\- `p0` (=2 by default): initial value of the factorization size p  
\- `AL_maxiter` (=20 by default): maximum number of augmented Lagrangian iterations   
\- `theta` (=1e-1 by default): threshold for estimating matrix ranks  
\- `delta` (=8 by default): number of negative eigenvalues used to construct a descent direction  
\- `alpha` (=0.5 by default): step size for escaping from saddle points  
\- `tolgradnorm` (=1e-8 by default): tolerance of stopping gradient norms used in the Riemannian trust-region method  
\- `TR_maxiter` (=40 by default): maximum number of Riemannian trust-region iterations  
\- `TR_maxinner` (=100 by default): maximum Hessian calls per trust-region iteration  
\- `line_search` (=0 by default): whether (=1) or not (=0) to use line search to decide the step size for escaping from saddle points  

**Output**  
`sol`: optimal solution  
`opt`: optimum  
`data`: auxiliary data

### SDPs with unit diagonal constraints
$$\inf_{X\succeq0}\langle C, X\rangle{\ }{\ }\text{s.t.}{\ }{\ }\mathcal{A}(X)=b,{\ }X_{ii}=1,{\ }\text{for}{\ }i=1,\ldots,n.$$

```matlab
clear options;
options.tol = 1e-8;
[sol, opt, data] = ManiSDP_unitdiag(At, b, c, K, options);
```

**Input**  
`At, b, c, K`: SeDuMi format data   
`options`:  
\- `tol` (=1e-8 by default): tolerance of maximal KKT residues  
\- `p0` (=2 by default): initial value of the factorization size p  
\- `AL_maxiter` (=300 by default): maximum number of augmented Lagrangian iterations  
\- `sigma0` (=1e-3 by default): initial value of the penalty parameter  
\- `sigma_min` (=1e-2 by default): minimum value of the penalty parameter  
\- `theta` (=1e-3 by default): threshold for estimating matrix ranks  
\- `delta` (=8 by default): number of negative eigenvalues used to construct a descent direction  
\- `alpha` (=0.1 by default): step size for escaping from saddle points  
\- `tolgradnorm` (=1e-8 by default): tolerance of stopping gradient norms used in the Riemannian trust-region method  
\- `TR_maxiter` (=4 by default): maximum number of Riemannian trust-region iterations  
\- `TR_maxinner` (=20 by default): maximum Hessian calls per trust-region iteration  
\- `tau1` (=1 by default): left factor for self-adaptively updating the penalty parameter  
\- `tau2` (=1 by default): right factor for self-adaptively updating the penalty parameter  
\- `line_search` (=0 by default): whether (=1) or not (=0) to use line search to decide the step size for escaping from saddle points  

**Output**  
`sol`: optimal solution  
`opt`: optimum  
`data`: auxiliary data

### SDPs with the unit trace constraint
$$\inf_{X\succeq0}\langle C, X\rangle{\ }{\ }\text{s.t.}{\ }{\ }\mathcal{A}(X)=b,{\ }\sum_{i}X_{ii}=1.$$

For the case that $X$ has a constant trace $c$, one can scale $X$ by the factor $\frac{1}{c}$ to match the unit trace case.

```matlab
clear options;
options.tol = 1e-8;
[sol, opt, data] = ManiSDP_unittrace(At, b, c, K, options);
```

**Input**  
`At, b, c, K`: SeDuMi format data   
`options`:  
\- `tol` (=1e-8 by default): tolerance of maximal KKT residues  
\- `p0` (=1 by default): initial value of the factorization size p  
\- `AL_maxiter` (=1000 by default): maximum number of augmented Lagrangian iterations  
\- `sigma0` (=1e1 by default): initial value of the penalty parameter  
\- `sigma_min` (=1e2 by default): minimum value of the penalty parameter  
\- `theta` (=1e-2 by default): threshold for estimating matrix ranks  
\- `delta` (=8 by default): number of negative eigenvalues used to construct a descent direction  
\- `alpha` (=0.05 by default): step size for escaping from saddle points  
\- `tolgradnorm` (=1e-8 by default): tolerance of stopping gradient norms used in the Riemannian trust-region method  
\- `TR_maxiter` (=3 by default): maximum number of Riemannian trust-region iterations  
\- `TR_maxinner` (=40 by default): maximum Hessian calls per trust-region iteration  
\- `tau1` (=1e-5 by default): left factor for self-adaptively updating the penalty parameter  
\- `tau2` (=1e-4 by default): right factor for self-adaptively updating the penalty parameter  
\- `line_search` (=1 by default): whether (=1) or not (=0) to use line search to decide the step size for escaping from saddle points  

**Output**  
`sol`: optimal solution  
`opt`: optimum  
`data`: auxiliary data

### SDPs with arbitrary affine constraints
$$\inf_{X\succeq0}\langle C, X\rangle{\ }{\ }\text{s.t.}{\ }{\ }\mathcal{A}(X)=b.$$

```matlab
clear options;
options.tol = 1e-8;
[sol, opt, data] = ManiSDP(At, b, c, K, options);
```

**Input**  
`At, b, c, K`: SeDuMi format data   
`options`:  
\- `tol` (=1e-8 by default): tolerance of maximal KKT residues  
\- `p0` (=1 by default): initial value of the factorization size p  
\- `AL_maxiter` (=1000 by default): maximum number of augmented Lagrangian iterations  
\- `sigma0` (=1e-2 by default): initial value of the penalty parameter  
\- `sigma_min` (=1e-1 by default): minimum value of the penalty parameter  
\- `theta` (=1e-2 by default): threshold for estimating matrix ranks  
\- `delta` (=8 by default): number of negative eigenvalues used to construct a descent direction  
\- `alpha` (=0.1 by default): step size for escaping from saddle points  
\- `tolgradnorm` (=1e-8 by default): tolerance of stopping gradient norms used in the Riemannian trust-region method  
\- `TR_maxiter` (=4 by default): maximum number of Riemannian trust-region iterations  
\- `TR_maxinner` (=20 by default): maximum Hessian calls per trust-region iteration  
\- `tau1` (=1e-2 by default): left factor for self-adaptively updating the penalty parameter  
\- `tau2` (=1e-1 by default): right factor for self-adaptively updating the penalty parameter  
\- `line_search` (=1 by default): whether (=1) or not (=0) to use line search to decide the step size for escaping from saddle points  

**Output**  
`sol`: optimal solution  
`opt`: optimum  
`data`: auxiliary data

### SDPs with multi-blocks
ManiSDP supports SDPs with multi-blocks:

$$\inf_{X}\langle C, X\rangle{\ }{\ }\text{s.t.}{\ }{\ }\mathcal{A}(X)=b, {\ }X \in \mathbb{S}_+^{n_1\times\cdots\times n_t}, {\ }\mathrm{diag}(X_i) = 1, {\ }i = 1,\ldots,K.nob.$$

```matlab
clear options;
options.tol = 1e-4;
K.nob = 10; % indicate that the first K.nob blocks are unit-diagonal
[sol, opt, data] = ManiSDP_multiblock(At, b, c, K, options);
```

**Input**  
`At, b, c, K`: SeDuMi format data   
`options`:  
\- `tol` (=1e-8 by default): tolerance of maximal KKT residues  
\- `AL_maxiter` (=1000 by default): maximum number of augmented Lagrangian iterations  
\- `sigma0` (=1e-2 by default): initial value of the penalty parameter  
\- `sigma_min` (=1e-1 by default): minimum value of the penalty parameter  
\- `theta` (=1e-2 by default): threshold for estimating matrix ranks  
\- `delta` (=8 by default): number of negative eigenvalues used to construct a descent direction  
\- `alpha` (=0.1 by default): step size for escaping from saddle points  
\- `tolgradnorm` (=1e-8 by default): tolerance of stopping gradient norms used in the Riemannian trust-region method  
\- `TR_maxiter` (=4 by default): maximum number of Riemannian trust-region iterations  
\- `TR_maxinner` (=20 by default): maximum Hessian calls per trust-region iteration  
\- `tau1` (=1e1 by default): left factor for self-adaptively updating the penalty parameter  
\- `tau2` (=1e1 by default): right factor for self-adaptively updating the penalty parameter  
\- `line_search` (=0 by default): whether (=1) or not (=0) to use line search to decide the step size for escaping from saddle points  

**Output**  
`sol`: optimal solution  
`opt`: optimum  
`data`: auxiliary data

### Solving second-order moment relaxations of binary quadratic programs
The second-order moment relaxation of the binary quadratic program

$$\inf_{x\in\mathbb{R}^q}x^{\intercal}Qx + e^{\intercal}x{\ }{\ }\mathrm{s.t.}{\ }{\ }x_i^2=1,i=1,\ldots,q$$

can be solved as follows.

```matlab
[At, b, c, K] = bqpmom(q, Q, e); % generate moment-SDP
options.tol = 1e-8;
[sol, opt, data] = ManiSDP_unitdiag(At, b, c, K, options);
```

### Solving moment relxations for polynomial optimization problems
$$\inf_{x\in\mathbb{R}^q}f(x){\ }{\ }\text{s.t.}{\ }{\ }h_i(x)=0,{\ }i=1,\ldots,l.$$

First, define your polynomial optimization problem using [SPOTLESS](https://github.com/spot-toolbox/spotless) as follows:
```matlab
q = 10;
x = msspoly('x', q);
f = x'*randn(q, 1);
h = x.^2 - 1;
problem.vars = x;
problem.objective = f;
problem.equality = h;
kappa = 2; % relaxation order
[SDP,info] = dense_sdp_relax(problem, kappa);
[sol, opt, data] = ManiSDP(At, b, c, n, options);
```

Then, generate the moment relxation and solve it with ManiSDP:
```matlab
kappa = 2; % relaxation order
[SDP,info] = dense_sdp_relax(problem, kappa);
options.tol = 1e-8;
[sol, opt, data] = ManiSDP(SDP.sedumi.At, SDP.sedumi.b, SDP.sedumi.c, SDP.sedumi.K, options);
```

## Solving low-rank SDPs with unit diagonal constraints in the dual form
ManiSDP also supports solving the following low-rank SDP using a dual Riemannian ADMM algorithm:

$$\inf_{y\in\mathbb{R}^m}b^{\intercal}y{\ }{\ }\text{s.t.}{\ }{\ }S=\mathcal{A}^{*}(y)+C\succeq0,\mathrm{diag}(S)=1,$$

where the linear operator $\mathcal{A}\mathcal{A}^*$ is assumed to be invertible. Here, low-rank means the SDP admits a low-rank optimal solution $S^{\star}$. In particular, SOS relaxations of binary quadratic programs are SDPs of the above form.

```matlab
clear options;
options.tol = 1e-8;
[sol, opt, data] = ManiDSDP_unitdiag(A, b, c, K, options);
```

**Input**  
`A, b, c, K`: SeDuMi format data   
`options`:  
\- `tol` (=1e-8 by default): tolerance of maximal KKT residues  
\- `p0` (=$\lceil log_2(m)\rceil$ by default): initial value of the factorization size p  
\- `ADMM_maxiter` (=300 by default): maximum number of ADMM iterations  
\- `sigma0` (=1e-3 by default): initial value of the penalty parameter  
\- `sigma_min` (=1e-3 by default): minimum value of the penalty parameter  
\- `theta` (=1e-3 by default): threshold for estimating matrix ranks  
\- `delta` (=8 by default): number of negative eigenvalues used to construct a descent direction  
\- `alpha` (=0.1 by default): step size for escaping from saddle points  
\- `tolgradnorm` (=1e-8 by default): tolerance of stopping gradient norms used in the Riemannian trust-region method  
\- `TR_maxiter` (=4 by default): maximum number of Riemannian trust-region iterations  
\- `TR_maxinner` (=20 by default): maximum Hessian calls per trust-region iteration  
\- `tau1` (=1e1 by default): left factor for self-adaptively updating the penalty parameter  
\- `tau2` (=1e2 by default): right factor for self-adaptively updating the penalty parameter  
\- `line_search` (=0 by default): whether (=1) or not (=0) to use line search to decide the step size for escaping from saddle points  

**Output**  
`sol`: optimal solution  
`opt`: optimum  
`data`: auxiliary data

### Solving second-order SOS relaxations of binary quadratic programs
The second-order SOS relaxation of the binary quadratic program
$$\inf_{x\in\mathbb{R}^q}x^{\intercal}Qx + e^{\intercal}x{\ }{\ }\mathrm{s.t.}{\ }{\ }x_i^2=1,i=1,\ldots,q$$
can be solved as follows.

```matlab
[A, b, dAAt, mb] = bqpsos(Q, e, q); % generate SOS-SDP
K.f = 1; 
K.s = mb;
c = [1; zeros(mb^2,1)];
v = zeros(size(A,1),1);
v(1) = 1;
A = [v A];
options.dAAt = dAAt;
options.tol = 1e-8;
options.line_search = 1;
[sol, opt, data] = ManiDSDP_unitdiag(A, b, c, K, options);
```

### Extension to the multi-block case
$$\inf_{y\in\mathbb{R}^m}b^{\intercal}y{\ }{\ }\text{s.t.}{\ }{\ }S=\mathcal{A}^*(y)+C\succeq0,\mathrm{diag}(S)=1,$$

```matlab
clear options;
options.tol = 1e-8;
K.nob = 10; % indicate that the first K.nob blocks are unit-diagonal
[sol, opt, data] = ManiDSDP_multiblock(A, b, c, K, options);
```

**Input**  
`A, b, c, K`: SeDuMi format data   
`options`:  
\- `tol` (=1e-8 by default): tolerance of maximal KKT residues  
\- `ADMM_maxiter` (=1000 by default): maximum number of ADMM iterations  
\- `sigma0` (=1e-1 by default): initial value of the penalty parameter  
\- `sigma_min` (=1e-2 by default): minimum value of the penalty parameter  
\- `theta` (=1e-2 by default): threshold for estimating matrix ranks  
\- `delta` (=8 by default): number of negative eigenvalues used to construct a descent direction  
\- `alpha` (=0.2 by default): step size for escaping from saddle points  
\- `tolgradnorm` (=1e-8 by default): tolerance of stopping gradient norms used in the Riemannian trust-region method  
\- `TR_maxiter` (=4 by default): maximum number of Riemannian trust-region iterations  
\- `TR_maxinner` (=20 by default): maximum Hessian calls per trust-region iteration  
\- `tau1` (=1e1 by default): left factor for self-adaptively updating the penalty parameter  
\- `tau2` (=1e1 by default): right factor for self-adaptively updating the penalty parameter  
\- `line_search` (=1 by default): whether (=1) or not (=0) to use line search to decide the step size for escaping from saddle points  

**Output**  
`sol`: optimal solution  
`opt`: optimum  
`data`: auxiliary data

## Examples
Check the folder `/example`.  

Note that to run `example_stls.m`, one has to first add the folder `NearestRankDeficient` in [CertifiablyRobustPerception](https://github.com/MIT-SPARK/CertifiablyRobustPerception); to run `example_rotationsearch.m`, one has to first add the folder `RotationSearch` in [CertifiablyRobustPerception](https://github.com/MIT-SPARK/CertifiablyRobustPerception).

## The Python version of ManiSDP
- [ManiSDP-Python](https://github.com/husisy/ManiSDP-matlab/tree/zc-dev)

## The Julia version of ManiSDP
Coming soon.

## References
[1] [Jie Wang and Liangbing Hu, Solving Low-Rank Semidefinite Programs via Manifold Optimization, 2025](https://link.springer.com/article/10.1007/s10915-025-02952-8)  
[2] [Jie Wang, Liangbing Hu, Bican Xia, A Dual Riemannian ADMM Algorithm for Low-Rank SDPs with Unit Diagonal, 2025](http://arxiv.org/abs/2512.04406)

## Contact
[Jie Wang](https://wangjie212.github.io/jiewang/): wangjie212@amss.ac.cn  
