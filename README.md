# ManiSDP
ManiSDP aims to solve the following low-rank semidefinite program (SDP) via manifold optimization:
$$\inf_{X\succeq0}{\ }\langle C, X\rangle{\ }\text{s.t.}{\ }\mathcal{A}(X)=b,{\ }\mathcal{B}(X)=d,$$
where the linear constraints $\mathcal{A}(X)=b$ are arbitrary while the linear constraints $\mathcal{B}(X)=d$ are assumed to impose certain manifold structures on the domain if present. Here, low-rank means the SDP admits a low-rank optimal solution.

## Dependencies
- [Manopt](https://github.com/NicolasBoumal/manopt)

## Usage
The optimal setting of parameters in ManiSDP is highly problem-dependent. We encourage the users to find the optimal parameters (defined in options) via preliminary experiments before running large-scale cases. So far, ManiSDP supports four types of SDPs.
### SDPs with only unit diagonal constraints
$$\inf_{X\succeq0}{\ }\langle C, X\rangle{\ }\text{s.t.}{\ }X_{ii}=1,{\ }\text{for}{\ }i=1,\ldots,n.$$

```matlab
clear options;
options.tol = 1e-8;
[sol, opt, data] = ManiSDP_onlyunitdiag(C, options);
```

Input  
C: cost matrix  
options:  
\- tol (=1e-8 by default): tolerance of residues
\- p0 (=2 by default): initial value of the factorization size p  
\- AL_maxiter (=20 by default): maximum number of iterations of the augmented Lagrangian method   
\- theta (=1e-1 by default): threshold for estimating matrix ranks  
\- delta (=8 by default): number of negative eigenvalues used to construct a descent direction  
\- alpha (=0.5 by default): step size for escaping from saddle points  
\- tolgrad (=1e-8 by default): tolerance of stopping gradient norms used in the Riemannian trust-region method  
\- TR_maxiter (=40 by default): maximum number of iterations of the Riemannian trust-region method  
\- TR_maxinner (=100 by default): maximum Hessian calls per trust-region iteration  
\- line_search (=0 by default): whether (=1) or not (=0) to use line search to decide the step size for escaping from saddle points  


Output  
sol: optimal solution  
opt: optimum  
data: auxiliary data

### SDPs with unit diagonal constraints
$$\inf_{X\succeq0}{\ }\langle C, X\rangle{\ }\text{s.t.}{\ }\mathcal{A}(X)=b,{\ }X_{ii}=1,{\ }\text{for}{\ }i=1,\ldots,n.$$

```matlab
clear options;
options.tol = 1e-8;
[sol, opt, data] = ManiSDP_unitdiag(At, b, c, n, options);
```

Input  
At, b, c: SeDuMi format data   
n: size of the PSD matrix  
options:  
\- tol (=1e-8 by default): tolerance of residues
\- p0 (=2 by default): initial value of the factorization size p  
\- AL_maxiter (=300 by default): maximum number of iterations of the augmented Lagrangian method  
\- sigma0 (=1e-3 by default): initial value of the penalty parameter  
\- sigma_min (=1e-2 by default): minimum value of the penalty parameter  
\- theta (=1e-3 by default): threshold for estimating matrix ranks  
\- delta (=8 by default): number of negative eigenvalues used to construct a descent direction  
\- alpha (=0.1 by default): step size for escaping from saddle points  
\- tolgrad (=1e-8 by default): tolerance of stopping gradient norms used in the Riemannian trust-region method  
\- TR_maxiter (=4 by default): maximum number of iterations of the Riemannian trust-region method  
\- TR_maxinner (=25 by default): maximum Hessian calls per trust-region iteration  
\- tao (=1 by default): factor for self-adaptively updating the penalty parameter  
\- line_search (=0 by default): whether (=1) or not (=0) to use line search to decide the step size for escaping from saddle points  

Output  
sol: optimal solution  
opt: optimum  
data: auxiliary data

### SDPs with the unit trace constraint
$$\inf_{X\succeq0}{\ }\langle C, X\rangle{\ }\text{s.t.}{\ }\mathcal{A}(X)=b,{\ }\sum_{i}X_{ii}=1.$$

For the case that $X$ has a constant trace $c$, one can scale $X$ by the factor $\frac{1}{c}$ to match the unit trace case.

```matlab
clear options;
options.tol = 1e-8;
[sol, opt, data] = ManiSDP_unittrace(At, b, c, n, options);
```

Input  
At, b, c: SeDuMi format data   
n: size of the PSD matrix  
options:  
\- tol (=1e-8 by default): tolerance of residues
\- p0 (=1 by default): initial value of the factorization size p  
\- AL_maxiter (=300 by default): maximum number of iterations of the augmented Lagrangian method  
\- sigma0 (=1e1 by default): initial value of the penalty parameter  
\- sigma_min (=1e2 by default): minimum value of the penalty parameter  
\- theta (=1e-2 by default): threshold for estimating matrix ranks  
\- delta (=10 by default): number of negative eigenvalues used to construct a descent direction  
\- alpha (=0.04 by default): step size for escaping from saddle points  
\- tolgrad (=1e-8 by default): tolerance of stopping gradient norms used in the Riemannian trust-region method  
\- TR_maxiter (=3 by default): maximum number of iterations of the Riemannian trust-region method  
\- TR_maxinner (=40 by default): maximum Hessian calls per trust-region iteration  
\- tao (=1/6e3 by default): factor for self-adaptively updating the penalty parameter  
\- line_search (=1 by default): whether (=1) or not (=0) to use line search to decide the step size for escaping from saddle points  

Output  
sol: optimal solution  
opt: optimum  
data: auxiliary data

### SDPs with arbitrary affine constraints
$$\inf_{X\succeq0}{\ }\langle C, X\rangle{\ }\text{s.t.}{\ }\mathcal{A}(X)=b.$$

```matlab
clear options;
options.tol = 1e-8;
[sol, opt, data] = ManiSDP(At, b, c, n, options);
```

Input  
At, b, c: SeDuMi format data   
n: size of the PSD matrix  
options:  
\- tol (=1e-8 by default): tolerance of residues
\- p0 (=1 by default): initial value of the factorization size p  
\- AL_maxiter (=300 by default): maximum number of iterations of the augmented Lagrangian method  
\- sigma0 (=1e-2 by default): initial value of the penalty parameter  
\- sigma_min (=1e-1 by default): minimum value of the penalty parameter  
\- theta (=1e-1 by default): threshold for estimating matrix ranks  
\- delta (=8 by default): number of negative eigenvalues used to construct a descent direction  
\- alpha (=0.2 by default): step size for escaping from saddle points  
\- tolgrad (=1e-8 by default): tolerance of stopping gradient norms used in the Riemannian trust-region method  
\- TR_maxiter (=4 by default): maximum number of iterations of the Riemannian trust-region method  
\- TR_maxinner (=50 by default): maximum Hessian calls per trust-region iteration  
\- tao (=0.25 by default): factor for self-adaptively updating the penalty parameter  
\- line_search (=1 by default): whether (=1) or not (=0) to use line search to decide the step size for escaping from saddle points  

Output  
sol: optimal solution  
opt: optimum  
data: auxiliary data

## References
[Jie Wang and Liangbing Hu, Solving Low-Rank Semidefinite Programs via Manifold Optimization, 2023]()  

## Contact
[Jie Wang](https://wangjie212.github.io/jiewang/): wangjie212@amss.ac.cn  
