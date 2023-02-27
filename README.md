# ManiSDP
ManiSDP aims to solve the following low-rank semidefinite program (SDP) via manifold optimization:
$$\inf_{X\succeq0}{\ }\langle C, X\rangle \text{s.t.} \mathcal{A}(X)=b,  \mathcal{B}(X)=d,$$
where the linear constraints $\mathcal{A}(X)=b$ are arbitrary while the linear constraints $\mathcal{B}(X)=d$ are assumed to impose certain manifold structures on the domain if present. Here, low-rank means the SDP admits a low-rank optimal solution.

## Dependencies
- [Manopt](https://github.com/NicolasBoumal/manopt)

## Usage
So far, ManiSDP supports
### SDPs with unit diagonal

```matlab
clear options;
options.tol = 1e-8;
[sol, obj, data] = ManiSDP_unitdiag(At, b, c, n, options);
```

Input:  
At, b, c: SeDuMi format data   
n: size of the PSD matrix  
options:

Output:  
sol: optimal solution  
obj: optimum  
data: auxiliary data

## References
[Jie Wang and Liangbing Hu, Solving Low-Rank Semidefinite Programs via Manifold Optimization, 2023]()  

## Contact
[Jie Wang](https://wangjie212.github.io/jiewang/): wangjie212@amss.ac.cn  
