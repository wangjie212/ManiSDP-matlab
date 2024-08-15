# ManiSDP
ManiSDP aims to solve the following low-rank semidefinite program (SDP) via manifold optimization:
$$\inf_{X\succeq0}{\ }\langle C, X\rangle{\ }\text{s.t.}{\ }\mathcal{A}(X)=b,{\ }\mathcal{B}(X)=d,$$
where the linear constraints $\mathcal{A}(X)=b$ are arbitrary while the linear constraints $\mathcal{B}(X)=d$, if present, are assumed to define certain manifold structure. Here, low-rank means the SDP admits a low-rank optimal solution.

## Dependencies
- [MATLAB](https://ww2.mathworks.cn/products/matlab.html?s_tid=hp_products_matlab)
- [Manopt](https://github.com/NicolasBoumal/manopt)

## Usage
**Note:** The optimal setting of parameters in ManiSDP is highly problem-dependent. We encourage the users to find the optimal parameters (defined in `options`) via preliminary experiments before running large-scale cases. In our experiences, the parameters (`sigma0`, `theta`, `TR_maxiter`, `TR_maxinner`, `tao`) have a significant influence on the performance of ManiSDP.  

ManiSDP accepts [SeDuMi](https://sedumi.ie.lehigh.edu/) format data.  

So far, ManiSDP supports four types of SDPs.

### SDPs with only unit diagonal constraints
$$\inf_{X\succeq0}{\ }\langle C, X\rangle{\ }\text{s.t.}{\ }X_{ii}=1,{\ }\text{for}{\ }i=1,\ldots,n.$$

```matlab
clear options;
options.tol = 1e-8;
[sol, opt, data] = ManiSDP_onlyunitdiag(C, options);
```

**Input**  
`C`: cost matrix  
`options`:  
\- `tol` (=1e-8 by default): tolerance of residues  
\- `p0` (=40 by default): initial value of the factorization size p  
\- `AL_maxiter` (=20 by default): maximum number of iterations of the augmented Lagrangian method   
\- `theta` (=1e-1 by default): threshold for estimating matrix ranks  
\- `delta` (=8 by default): number of negative eigenvalues used to construct a descent direction  
\- `alpha` (=0.5 by default): step size for escaping from saddle points  
\- `tolgrad` (=1e-8 by default): tolerance of stopping gradient norms used in the Riemannian trust-region method  
\- `TR_maxiter` (=40 by default): maximum number of iterations of the Riemannian trust-region method  
\- `TR_maxinner` (=100 by default): maximum Hessian calls per trust-region iteration  
\- `line_search` (=0 by default): whether (=1) or not (=0) to use line search to decide the step size for escaping from saddle points  

**Output**  
`sol`: optimal solution  
`opt`: optimum  
`data`: auxiliary data

### SDPs with unit diagonal constraints
$$\inf_{X\succeq0}{\ }\langle C, X\rangle{\ }\text{s.t.}{\ }\mathcal{A}(X)=b,{\ }X_{ii}=1,{\ }\text{for}{\ }i=1,\ldots,n.$$

```matlab
clear options;
options.tol = 1e-8;
[sol, opt, data] = ManiSDP_unitdiag(At, b, c, n, options);
```

**Input**  
`At, b, c`: SeDuMi format data   
`n`: size of the PSD matrix  
`options`:  
\- `tol` (=1e-8 by default): tolerance of residues  
\- `p0` (=2 by default): initial value of the factorization size p  
\- `AL_maxiter` (=300 by default): maximum number of iterations of the augmented Lagrangian method  
\- `sigma0` (=1e-3 by default): initial value of the penalty parameter  
\- `sigma_min` (=1e-2 by default): minimum value of the penalty parameter  
\- `theta` (=1e-3 by default): threshold for estimating matrix ranks  
\- `delta` (=8 by default): number of negative eigenvalues used to construct a descent direction  
\- `alpha` (=0.1 by default): step size for escaping from saddle points  
\- `tolgrad` (=1e-8 by default): tolerance of stopping gradient norms used in the Riemannian trust-region method  
\- `TR_maxiter` (=4 by default): maximum number of iterations of the Riemannian trust-region method  
\- `TR_maxinner` (=25 by default): maximum Hessian calls per trust-region iteration  
\- `tao` (=1 by default): factor for self-adaptively updating the penalty parameter  
\- `line_search` (=0 by default): whether (=1) or not (=0) to use line search to decide the step size for escaping from saddle points  

**Output**  
`sol`: optimal solution  
`opt`: optimum  
`data`: auxiliary data

### SDPs with the unit trace constraint
$$\inf_{X\succeq0}{\ }\langle C, X\rangle{\ }\text{s.t.}{\ }\mathcal{A}(X)=b,{\ }\sum_{i}X_{ii}=1.$$

For the case that $X$ has a constant trace $c$, one can scale $X$ by the factor $\frac{1}{c}$ to match the unit trace case.

```matlab
clear options;
options.tol = 1e-8;
[sol, opt, data] = ManiSDP_unittrace(At, b, c, n, options);
```

**Input**  
`At, b, c`: SeDuMi format data   
`n`: size of the PSD matrix  
`options`:  
\- `tol` (=1e-8 by default): tolerance of residues  
\- `p0` (=1 by default): initial value of the factorization size p  
\- `AL_maxiter` (=300 by default): maximum number of iterations of the augmented Lagrangian method  
\- `sigma0` (=1e1 by default): initial value of the penalty parameter  
\- `sigma_min` (=1e2 by default): minimum value of the penalty parameter  
\- `theta` (=1e-2 by default): threshold for estimating matrix ranks  
\- `delta` (=10 by default): number of negative eigenvalues used to construct a descent direction  
\- `alpha` (=0.04 by default): step size for escaping from saddle points  
\- `tolgrad` (=1e-8 by default): tolerance of stopping gradient norms used in the Riemannian trust-region method  
\- `TR_maxiter` (=3 by default): maximum number of iterations of the Riemannian trust-region method  
\- `TR_maxinner` (=40 by default): maximum Hessian calls per trust-region iteration  
\- `tao` (=1/6e3 by default): factor for self-adaptively updating the penalty parameter  
\- `line_search` (=1 by default): whether (=1) or not (=0) to use line search to decide the step size for escaping from saddle points  

**Output**  
`sol`: optimal solution  
`opt`: optimum  
`data`: auxiliary data

### SDPs with arbitrary affine constraints
$$\inf_{X\succeq0}{\ }\langle C, X\rangle{\ }\text{s.t.}{\ }\mathcal{A}(X)=b.$$

```matlab
clear options;
options.tol = 1e-8;
[sol, opt, data] = ManiSDP(At, b, c, n, options);
```

**Input**  
`At, b, c`: SeDuMi format data   
`n`: size of the PSD matrix  
`options`:  
\- `tol` (=1e-8 by default): tolerance of residues  
\- `p0` (=1 by default): initial value of the factorization size p  
\- `AL_maxiter` (=300 by default): maximum number of iterations of the augmented Lagrangian method  
\- `sigma0` (=1e-2 by default): initial value of the penalty parameter  
\- `sigma_min` (=1e-1 by default): minimum value of the penalty parameter  
\- `theta` (=1e-1 by default): threshold for estimating matrix ranks  
\- `delta` (=8 by default): number of negative eigenvalues used to construct a descent direction  
\- `alpha` (=0.2 by default): step size for escaping from saddle points  
\- `tolgrad` (=1e-8 by default): tolerance of stopping gradient norms used in the Riemannian trust-region method  
\- `TR_maxiter` (=4 by default): maximum number of iterations of the Riemannian trust-region method  
\- `TR_maxinner` (=50 by default): maximum Hessian calls per trust-region iteration  
\- `tao` (=0.25 by default): factor for self-adaptively updating the penalty parameter  
\- `line_search` (=1 by default): whether (=1) or not (=0) to use line search to decide the step size for escaping from saddle points  

**Output**  
`sol`: optimal solution  
`opt`: optimum  
`data`: auxiliary data

### SDPs with multi-blocks
ManiSDP supports SDPs with multi-blocks:

$$\inf_{X}{\ }\langle C, X\rangle{\ }\text{s.t.}{\ }\mathcal{A}(X)=b, {\ }X \in \mathbb{S}_+^{n_1\times\cdots\times n_t}, {\ }\mathrm{diag}(X_i) = 1, {\ }i = 1,\ldots,K.nob.$$

```matlab
clear options;
options.tol = 1e-4;
K.nob = 10; % indicate the first K.nob blocks are unit-diagonal
[sol, opt, data] = ManiSDP_multiblock(At, b, c, K, options);
```

### Solving moment relxations for polynomial optimization problems
$$\inf_{x\in\mathbb{R}^q}{\ }f(x){\ }\text{s.t.}{\ }h_i(x)=0,{\ }i=1,\ldots,l.$$

First define your polynomial optimization problem using [SPOTLESS](https://github.com/spot-toolbox/spotless) as follows:
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

Then generate the moment relxation and solve it with ManiSDP:
```matlab
kappa = 2; % relaxation order
[SDP,info] = dense_sdp_relax(problem, kappa);
options.tol = 1e-8;
[sol, opt, data] = ManiSDP(SDP.sedumi.At, SDP.sedumi.b, SDP.sedumi.c, SDP.sedumi.K.s, options);
```

## Examples
Check the folder `/example`.  

Note that to run `example_stls.m`, one has to first add the folder `NearestRankDeficient` in [CertifiablyRobustPerception](https://github.com/MIT-SPARK/CertifiablyRobustPerception); to run `example_rotationsearch.m`, one has to first add the folder `RotationSearch` in [CertifiablyRobustPerception](https://github.com/MIT-SPARK/CertifiablyRobustPerception).

## The Julia version for ManiSDP
Coming soon.

## References
[Jie Wang and Liangbing Hu, Solving Low-Rank Semidefinite Programs via Manifold Optimization, 2023](http://arxiv.org/abs/2303.01722)  

## Contact
[Jie Wang](https://wangjie212.github.io/jiewang/): wangjie212@amss.ac.cn  
