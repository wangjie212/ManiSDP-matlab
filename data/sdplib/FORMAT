The SDP problems in this directory are all encoded in the SDPA sparse
format [1].  This file provides a description of the file format.  

The SDP Problem:

We work with a semidefinite programming problem that has been written
in the following standard form:
 
(P)    min c1*x1+c2*x2+...+cm*xm
       st  F1*x1+F2*x2+...+Fm*xn-F0=X
                                 X >= 0

The dual of the problem is:
 
(D)    max tr(F0*Y)
       st  tr(Fi*Y)=ci           i=1,2,...,m
                 Y >= 0

Here all of the matrices F0, F1, ..., Fm, X, and Y are assumed to be
symmetric of size n by n.  The constraints X>=0 and Y>=0 mean that X
and Y must be positive semidefinite.  

Note that several other standard forms for SDP have been used by 
a number of authors- these can generally be translated into the SDPA
standard form with little effort.  


File Format:

The file consists of six sections:
 
1. Comments.  The file can begin with arbitrarily many lines of comments.
Each line of comments must begin with '"' or '*'.  

2. The first line after the comments contains m, the number of constraint
matrices.  Additional text on this line after m is ignored.
 
3. The second line after the comments contains nblocks, the number of 
blocks in the block diagonal structure of the matrices.  Additional text
on this line after nblocks is ignored.  
 
4. The third line after the comments contains a vector of numbers that 
give the sizes of the individual blocks.  The special characters 
',', '(', ')', '{', and '}' can be used as punctuation and are ignored.  
Negative numbers may be used to indicate that a block is actually a diagonal
submatrix.  Thus a block size of "-5" indicates a 5 by 5 block in which 
only the diagonal elements are nonzero.  

5. The fourth line after the comments contains the objective function
vector c.  
 
6. The remaining lines of the file contain entries in the constraint
matrices, with one entry per line.  The format for each line is 
 
  <matno> <blkno> <i> <j> <entry>
 
Here <matno> is the number of the matrix to which this entry belongs, 
<blkno> specifies the block within this matrix, <i> and <j> specify a
location within the block, and <entry> gives the value of the entry in
the matrix.  Note that since all matrices are assumed to be symmetric, 
only entries in the upper triangle of a matrix are given.  

Example:

         min 10x1+20x2

         st  X=F1x1+F2x2-F0
 
             X >= 0
 
where
 
 
  F0=[1 0 0 0
      0 2 0 0
      0 0 3 0
      0 0 0 4]

  F1=[1 0 0 0
      0 1 0 0
      0 0 0 0
      0 0 0 0]

  F2=[0 0 0 0
      0 1 0 0
      0 0 5 2
      0 0 2 6]


In SDPA sparse format, this problem can be written as:
 
"A sample problem.  
2 =mdim
2 =nblocks
{2, 2}
10.0 20.0
0 1 1 1 1.0
0 1 2 2 2.0
0 2 1 1 3.0
0 2 2 2 4.0
1 1 1 1 1.0
1 1 2 2 1.0
2 1 2 2 1.0
2 2 1 1 5.0
2 2 1 2 2.0
2 2 2 2 6.0

 
References:
 
[1] K. Fujisawa, M. Kojima, and K. Nakata.  SDPA (Semidefinite Programming 
Algorithm) User's Manual.  Technical Report B-308, Department of Mathematical
and Computing Sciences, Tokyo Institute of Technology. Revised, May, 1998.
