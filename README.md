

# <span style="color:blue">fastLinearAlgebra</span>

Quick REVIEW of linear algebra. 

## Author:

Norm Matloff <br>
Prof. of Computer Science <br>
UC Davis <br>
[my bio](http://heather.cs.ucdavis.edu/matloff.html)

# Who Is This for?

The main target audience is data science students and professionals who
have some previous background in linear algebra but now need a review.

The R language is used for examples and discussion on use of R in linear
algebra, but is not necessary for following the mathematical
presentation.

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

You probably studied *simultaneous linear equations* in high school,
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

in which the left-hand side consists of multiplication of a vector v by
a matrix M.  Here

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
w, if considered as matrices, are each of size of 2x1.

Note: We are displaying the vectors v and w here as columns, but
sometimes they will be shown as rows.  Either way, the key point is that
a vector is an ordered set of numbers, however it is displayed.

Now, what does that matrix equation mean?  First, consider the
expression Mv.  Motivated by our comments above regarding the
correspondence between linear equations and inner prodcts, we define a
matrix-vactor product this way:

> Mv is a new vector, of the same length as v.  The i<sup>th</sup>
> number in the new vector is the inner product with row i of M with v.

In R:

``` r
mTimesv <- function(M,v) 
{
   if (ncol(M) != length(v)) stop('incompatible sizes')
   # prepare space for the output vector
   mvProduct <- vector(length=length(v))
   # build the output vector
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

So, the matrix-vector product Mv for M and v as above is equal to

$$
\left (
\begin{array}{r}
3a+2b \\
2a - 2b \\
\end{array}
\right )
$$

and our equation Mv = w is

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
Now we will in turn build matrix-matrix multiplication by using
matrix-vector multiplication as the base.

> The matrix-matrix product MQ is a new matrix having the same number of
> rows as M and the same number of columns as Q.  The j<sup>th</sup> 
> column of a matrix-matrix product MQ is the matrix-vector product
> of M with the j<sup>th</sup> column of Q.

In R:

``` r

# matrix times matrix
mTimesq <- function(M,Q) 
{
   if (ncol(M) != nrow(Q)) stop('nonconformable factors')
   # prepare space for building the output matrix
   mq <- matrix(nrow=nrow(M),ncol=ncol(Q))
   # build output matrix
   for (j in 1:ncol(Q)) {
      mq[,j] <- mTimesv(M,Q[,j])
   }
   mq
}

```

## Or, build from vector-times-matrix

By the way, alternatively we could have first defined wR, i.e.
vector-matrix multiplication, yielding another vector whose
j<sup>th</sup> element is the dot product of w and column j of R.
We could then build on this to define matrix-matrix multiplication LR.

Either way, it's very important that the reader view the process in this
way, called *partitioned matrices*.  

## Partitioned matrices: an invaluable visualization tool 

Look at our matrix-vector product Mv.  Let's denote row i of M by
M<sub>i.</sub>.  Then our building of the matrix-vector product 
can be displayed--and most important, visualized--symbolicly. 
For an r-row matrix M and a conformable vector v, we think of Mv as

$$
\left (
\begin{array}{r}
M_{1.} v \\
M_{2.}  v\\
... \\
M_{r.}  v\\
\end{array}
\right )
$$

where M<sub>i.</sub>v is interpreted as the dot product of
M<sub>i.</sub> with v.

The expression
$$
\left (
\begin{array}{r}
M_{1.} \\
M_{2.} \\
... \\
M_{r.} \\
\end{array}
\right ) 
$$

is a column of r things, so if you squint, it "looks" like a column of r
numbers.  It wil turn out later that this way to visualize things is
useful.


Similarly, for the matrix-matrix product MQ, our above formulation can
be written symbolicly as follows, with Q<sub>.j</sub> symbolizing column
j of Q, j = 1,...,c, we think of MQ as 

$$
\left (
M Q_{.1},...,M Q_{.c}
\right )
$$

Mastering this way of visualizing matrix projects will go a long way to
helping you master linear algebra.

## "Real" matrix multiplication in R, etc.

Though the function mtimesq() above does produce the correct answer, it
is not efficient.  Instead, we use a special operator, %*%.

### %*%

``` r
 A <- rbind(4:5,8:9)
> A
     [,1] [,2]
[1,]    4    5
[2,]    8    9
> B <- rbind(c(1.5,0),c(1,1))
> B
     [,1] [,2]
[1,]  1.5    0
[2,]  1.0    1
> A %*% B
     [,1] [,2]
[1,]   11    5
[2,]   21    9
```

The %*% operator is written in C, so it is faster.  So, use it rather
than mtimesq(), but keep the code structuroe of the latter in mind for
visualization purposes.

### Internal storage

A computer's memory is linear, i.e. Word 0, then Word 1, Word 2 and so
on.  (If our data consists of numbers, one number is stored per word.
For character data, m characters are stored per word, where m is the
word size in bytes.)  But a matrix is two-dimensional, i.e. nonlinear.
So, how is a matrix stored?

The two main storage schemes are *column-major order* and 
*row-major order*.  C uses the latter, while R uses the former,
meaning that first all of column 1 is stored, then column 2, etc.

So for instance the matrix A above is stored as 4, 8, 5, 9.  Indeed, we
could have created A using that knowledge:

``` r
> matrix(c(4,8,5,9),ncol=2)
     [,1] [,2]
[1,]    4    5
[2,]    8    9
```

### Vectors as matrices, and vice versa

Consider this multiplication:

``` r
 A %*% c(3,2)
     [,1]
[1,]   22
[2,]   42
```

Even though the (3,2) was not a matrix, R treated it as a 2x1 matrix
(it would not be conformable if it were treated as 1x2),
and it did the multiplication.  And look at this one:

``` r
> c(3,2) %*% A
     [,1] [,2]
[1,]   28   33
```

R treated the (3,2) as a 1x2 matrix.

And look at this:

``` r
> A[3]
[1] 5
```

Here R treated A, which is 2x2, as a 4-element vector.

# Lesson 3: Linear Combinations of Vectors, Matrix Rank

In data science, a common approach to dealing with very large datasets
is to find a *low-rank approximation*.  Say we have medical data, with
450 variables for each patient (height, weight, age, waking pulse rate,
disease history etc.).  Analysis with so many variables might be
unweildy.  We cannot go into details here, but matrices are involved,
and one type would be 450x450.  A low-rank version, say rank m, and
would reduce the size of the matrix to mxm.  So, what is rank?

First, we need the notion of a *linear combination* of vectors, which is
a sum of scalar products of vectors.  E.g.

3 (3,2) - 0.5 (12,8) = (3,-1)

The 3 and -0.5 are the *coefficients* in the linear combination.

That leads to the notion of *matrix rank*.  Consider the matrix M:

$$
M = \left (
\begin{array}{rrr}
3 & 2 & 0.5 \\
2 & -2 & 0 \\
4 & 1 & 0.5 \\
\end{array}
\right )
$$

Row 3 of M happens to be 1 times Row 1 + 0.5 times Row 2.  So, Row 3 is
a linear combination of Rows 1 and 2.  And put another way,

1 Row 1 + 0.5 Row 2 - Row 3 = (0,0,0)

We say that the rows of M are *linearly dependent*, meaning that there
exists some linear combination of its rows that equals the 0 vector
(vector consisting of all 0s), where at least one coefficient of the
linaar combination is nonzero.  If no such linear combination were to
exist, we would say the rows are *linearly independent*.

In this case, although the set of three rows is not linearly
independent, any two of them are so.  We say that the rank of this
matrix is 2.

> The rank of a matrix (not necessarily a square matrix) is the maximal number
> of linearly independent rows.

Actually it turns out that we can define in terms of linearly
independent columns too--the row rank and column rank of a matrix can be
proven to be equal.  One implication of that is:

> The rank of an n x m matrix is at most min(n,m).

A matrix is said to be of *full rank* if its rank achieves that upper
bound.  E.g. if a matrix is of size 6x4 and has rank 4, it is of full
rank.

Another important inequality is

> The rank of the matrix product AB is at most the minimum of the rank
> of A and the rank of B.

# Lesson 4: Matrix Transpose and Inverse

## Transpose

> Given any (not necessary square) matrix A, its *transpose* is the matrix
> obtained by exchanging rows and columns.  We use ' or t to denote
> transpose, i.e. A' or A<sup>t</sup>.

In R, the t() function implements transpose.

``` r
> W <- matrix(1:6,nrow=2)
> W
     [,1] [,2] [,3]
[1,]    1    3    5
[2,]    2    4    6
> t(W)
     [,1] [,2]
[1,]    1    2
[2,]    3    4
[3,]    5    6
```

## Diagonal matrices, identity matrix

A useful function in R is diag(), which produces a square matrix with 0s
in the off-diagonal positions:

``` r
> diag(c(1,2,5))
     [,1] [,2] [,3]
[1,]    1    0    0
[2,]    0    2    0
[3,]    0    0    5
```

An important special case is that of the *identity matrix*, generally
denoted by I, which acts like a multiplicative 1 in numbers:

``` r
> I <- diag(3)
> I
     [,1] [,2] [,3]
[1,]    1    0    0
[2,]    0    1    0
[3,]    0    0    1
> B <- matrix(c(3,1,0,0.1,1,1,2:4),ncol=3)
> B
     [,1] [,2] [,3]
[1,]    3  0.1    2
[2,]    1  1.0    3
[3,]    0  1.0    4
> I %*% B
     [,1] [,2] [,3]
[1,]    3  0.1    2
[2,]    1  1.0    3
[3,]    0  1.0    4
> B %*% I
     [,1] [,2] [,3]
[1,]    3  0.1    2
[2,]    1  1.0    3
[3,]    0  1.0    4
```

Even though we speak of "the" identity matrix, of course there are
different ones for each matrix size, i.e. 2x2, 3x3 etc.

## Inverse

Since matrix algebra has something analogous to the number 1, maybe it
also has something like "division"?  For some matrices, yes.

> For a given square matrix A, if there is another square matrix B such
> that AB = I (which actually will also imply that BA = I), we say that B
> is the *inverse* of A, and write B = A<sup>-1</sup>.  And, A will have
> an inverse if and only if it is of full rank.

So A<sup>-1</sup> is roughly analogous to "1/A" in the numbers world.
Square matrices of less than full rank have no inverse, just as 0 has no
reciprocal in the numbers world.  Such matrices are called *singular*.

In R, we can obtain the inverse of a matrix using the solve() function:

``` r
> C <- rbind(1:3,c(2,0,0),c(1,5,12))
> C
     [,1] [,2] [,3]
[1,]    1    2    3
[2,]    2    0    0
[3,]    1    5   12
> D <- solve(C)
> D
           [,1]       [,2]       [,3]
[1,]  0.0000000  0.5000000  0.0000000
[2,]  1.3333333 -0.5000000 -0.3333333
[3,] -0.5555556  0.1666667  0.2222222
> C %*% D
     [,1]         [,2]         [,3]
[1,]    1 1.110223e-16 1.110223e-16
[2,]    0 1.000000e+00 0.000000e+00
[3,]    0 0.000000e+00 1.000000e+00
```

Note that we got a bit of roundoff error.  We can't avoid this entirely,
but there are advanced functions that reduce it.

In our old formulation of a set of linear equations as Mv = w,
we see now how to solve for v:

v = M<sup>-1</sup> w

# Lesson 5: Eigenvalues and Eigenvectors

> Let A be a square matrix.  If there is a nonzero vector v and a number 
> &lambda; such that Av = &lambda; v, we say that &lambda; is an 
> *eigenvalue* of A, with *eigenvector* v.

Note: "Nonzero vector means "not all elements are 0s.  The zero vector
is one in which all elements are 0s.

This notion also plays a key role in the low-rank approximation
processes we mentioned earlier as being key to data science.

We can use the eigen() function for this, e.g.:

``` r
> eigen(C)
eigen() decomposition
$values
[1] 12.4839914  1.4861788 -0.9701703

$vectors
          [,1]       [,2]       [,3]
[1,] 0.2592994 -0.5462472  0.4164962
[2,] 0.0415411 -0.7351030 -0.8586044
[3,] 0.9649032  0.4015440  0.2988801
> C %*% eigs$vectors[,1]
           [,1]
[1,]  3.2370912
[2,]  0.5185987
[3,] 12.0458432
> eigs$values[1] * eigs$vectors[,1]
[1]  3.2370912  0.5185987 12.0458432
# yes, they are equal
```

# Lessons 6 and Beyond

(Under construction.)

