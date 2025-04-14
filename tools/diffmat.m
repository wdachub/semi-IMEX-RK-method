function D=diffmat(x,n_s,m)
% this program is design in order to get the discrete derivative operator by finite difference methods with n_s points
% stencil. We use the inner point to get derivative at boundary point, do
% not use the derivative boundary condition.

%input:
%     x: the coordinate of stencil
%     n_s: number of stencil
%     m: the desired order of derivative. Requires length(x) > m & n_s>m
%oupt:
%     D: a cell which record all requied differentiation  Matrix.
%        D{1} First Derivative Matrix
%        D{n} n-th Derivative Matrix
n=length(x);
ns2=floor(n_s/2);

if m >= n
    error('length(x) must be larger than m')
elseif m > n_s
    error('n_s must be larger than m')
end
if n_s > n
    error('number of grid points must be larger than number of stencil points')
end

for i=m:-1:1
D{i} = sparse(n, n);
end

for i=1:n
    if (i <ns2+1)
        ii=(1:n_s);%The array containing the indices of stencil points for x(i)
        x2=x(ii);%The coordinate of the stencil points
    elseif (i>n-ns2)
        ii=n-n_s+1:n;%The array containing the indices of stencil points for x(i)
        x2=x(ii);%The coordinate of the stencil points
    else
        if mod(n_s,2)==1% n_s is odd.
            ii=(i-ns2:i+ns2);%The array containing the indices of stencil points for x(i)
            x2=x(ii);%The coordinate of the stencil points
        else
            ii=(i-ns2+1:i+ns2);%The array containing the indices of stencil points for x(i)
            x2=x(ii);%The coordinate of the stencil points
        end
    end
    c=fdcoeffF(x(i),x2,n_s,m);% m*n_s matrix containing the weights of of each stencil point
    for om=m:-1:1
        D{om}(i,ii)=c(om+1,1:n_s);%Elements of the Difrentiation matrix for row i and columbs ii, First Derivative
    end
end

end

function c=fdcoeffF(xbar,x_s,n_s,m)

% Compute coefficients for finite difference approximation for the
% derivative of order 0 to k at xbar based on grid values at points in x.

% OUTPUT:
% This function returns a matrix c of dimension k by n, where n=length(x),
% containing coefficients to approximate u^{(k)}(xbar),
% the zeroth to k'th derivative of u evaluated at xbar,  based on n_s values
% of u at x(1), x(2), ... x(n).
%
% If U is a column vector containing u(x) at these n points, then
% c(k,:)*U will give the approximation to u^{(k)}(xi).
%
% Note for k=0 this can be used to evaluate the interpolating polynomial
% itself.
%
%
%INPUT:
% m£º the desired order of derivative. Requires length(x) > m & n_s>m.
% xbar£º location at which we want to approximate the derivative (may but
% need not be a grid point).
% x£º the stencil points. row vector or column vector.
% Usually the elements x(i) are monotonically increasing
% and x(1) <= xbar <= x(n), but neither condition is required.
% The x values need not be equally spaced but must be distinct.
%n_s: the number of stencil points.
%
% This program should give the same results as fdcoeffV.m, but for large
% values of n is much more stable numerically.
%
% Based on the program "weights" in
%   B. Fornberg, "Calculation of weights in finite difference formulas",
%   SIAM Review 40 (1998), pp. 685-691.
%
% Note: Forberg's algorithm can be used to simultaneously compute the
% coefficients for derivatives of order 0, 1, ..., m where m <= n-1.
% This gives a coefficient matrix C(1:n,1:m) whose k'th column gives
% the coefficients for the k'th derivative.

% The following are some notes regarding the implementation of the algorithm:
% - Although the algorithm can be shown to be numerically robust, applying a
% difference stencil to smooth data can lead to cancellation of significant digits
% (especially when approximating high derivatives).
% - A call to weights to obtain the coefficients for the mth derivative returns also
% the coefficients for the kth derivative, k = 0, 1, ... ,m. These are "byproducts,"
% available at no additional computational cost.
% - The code can be modified to return all the data above also for stencils which
% extend only over xo, xl,... , xi, j = 0, 1,... , n-still at no additional cost
% [4]. (This version is actually slightly shorter than the present one; c is then a
% 3-D array, and less care needs to be taken to prevent premature overwritings
% within it.)
% - Calling weights with m = 0 gives the weights for polynomial interpolation
% at a cost (to leading order) of 2n2 floating point operations. This can be
% compared to 1.5n^2 operations to obtain Newton's divided differences and
% 2.5nr2 operations for interpolation using the well-known algorithms by Aitken
% and Neville. Particularly large savings are realized if several functions are to
% be interpolated on the same grid; the weights can then be reused at a cost of
% only 2n operations per case (this situation will arise already for a single data
% function if it is given on a 2-D or 3-D Cartesian grid).
%
n=length(x_s);
if m >=n
    error('length(x) must be larger than m')
elseif m > n_s
    error('n_s must be larger than m')
end
if n_s > n
    error('number of grid points must be larger than number of stencil points')
end
m=m+1; %including the zero-th derivative
c=zeros(n_s,m); c(1,1)=1.0;
c1=1.0;
c4=x_s(1)-xbar;

for i=2:n_s
    mn=min(i,m);
    c2=1.0;
    c5=c4;
    c4=x_s(i)-xbar;
    
    for j=1:i-1
        c3=x_s(i)-x_s(j);
        c2=c3*c2;
        for k=mn:-1:2
            c(i,k)=c1*((k-1)*c(i-1,k-1)-c5*c(i-1,k))/c2;
        end
        c(i,1)=-c1*c5*c(i-1,1)/c2;
        for k=mn:-1:2
            c(j,k)=(c4*c(j,k)-(k-1)*c(j,k-1))/c3;
        end
        c(j,1)=c4*c(j,1)/c3;
    end
    c1=c2;
end
c=c';
end