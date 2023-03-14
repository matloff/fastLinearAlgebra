

# <span style="color:blue">fastLinearAlgebra</span>

Quick review of linear algebra. Some facility with R helpful but not required.

### Norm Matloff, Prof. of Computer Science, UC Davis; [my bio](http://heather.cs.ucdavis.edu/matloff.html)

# Who Is This for?

The main target audience is data science students and professionals who
have some previous background in linear algebra but now need a review.

The R language is used for examples but is not necessary for following
the presentation.

# Lesson 1: Vectors and Dot Products

For now, we define a *vector* as an ordered set of numbers,
such (5,2.2) and (-1.12,6,88,8.8).  Note that this definition holds for
R as well, except that in R the components can be non-numeric.  That
latter situation won't concern us here, but it will always be important
to distinguish between a math definition and a coding one.

The *dot product* of two vectors (which must be of the same length) is
the sum of their elementwise products.  

Example:   The dot product of (2,5,1.2) and (0,3.3,8) is

2 x 0 + 5 x 3.3 + 1.2 x 8 = 26.1

In R:

``` r
dotProd <- function(x,y) {
   if (length(x) != length(y)) stop('x,y must be of the same length')
   sum(x*y)
}

# example
dotProd(c(3,2),c(-1.5,4))
# 3.5

```

# Lesson 2: Matrices and Matrix Products

## Matrix-vector products

You probably studies *simultaneous linear equations* in high school,
e.g.

3a + 2b = 8
<br>
a - 2b = 4

But now the point is that, say, that first equation can be written as

> dot product of (3,2) with (a,b) = 8

and the second as 

> dot product of (1,-2) with (a,b) = 4

This motivates the *matrix equation*

M v = w

where

$$
M = \left (
\begin{array}{rr}
3 & 2 \\
2 & -2 \\
\end{array}
\right )
$$

$$
v =\left (
\begin{array}{r}
a \\
b \\
\end{array}
\right )
$$

$$
w =\left (
\begin{array}{r}
8 \\
4 \\
\end{array}
\right )
$$

The matrix M is said to be of size 2x2 (2 rows, 2 columns), while v and
w are each matrices of size of 2x2.

Now, what does that matrix equation mean?  First, consider the
expression Mv.  It is defined this way:

> The number in the i<sup>th</sup> row of Mv is the inner product with row i of
> M with v.

In R:

``` r
mTimesv <- function(M,v) 
{
   if (ncol(M) != length(v)) stop('incompatible sizes')
   mvProduct <- vector(length=length(v))
   for (i in 1:nrow(M)) 
      mvProduct[i] <- dotProd(M[i,],v)
   mvProduct
}
#example
u <- rbind(c(1,5),c(2,2))
u
#     1    5
#     2    2
w <- c(5.5,3)
mTimesv(u,v)
# 20.5 17.0
```

So, the matrix-vector product Mv is equal to

$$
\left (
\begin{array}{r}
3a+2b \\
2a - 2b \\
\end{array}
\right )
$$

So, our equation Mv = w is

$$
\left (
\begin{array}{r}
3a+2b \\
2a - 2b \\
\end{array}
\right ) = 
\left (
\begin{array}{r}
8 \\
4 \\
\end{array}
\right ) 
$$

Setting equal the corresponding components of each side, we have

3a + 2b = 8
<br>
a - 2b = 4

i.e. our original set of equations.  So, the matrix forumulation Mv = w
is indeed equivalent to the original equations.



**The reader should make absolutely sure that he/she understands (a) how in
this manner, matrices provide us with a way to compactly describe
systems of linear equations, and (b) the role of dot products in
defining matrix-vector multiplication.**

Simple, yes, but really, no kidding; don't go further until you
understand things well--meaning you could explain it to others.

## Matrix-matrix products

We built matrix-vector multiplication by using dot product as the base.
Now we will build matrix-matrix multiplication by using matrix-vector
multiplication as the base.

> The j<sup>th</sup> column of a matrix product MQ is the product of M
> with the j<sup>th</sup> column of Q.

