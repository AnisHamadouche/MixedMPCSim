function y = sumsqr( x, dim )

%sumsqr   Sum of squares.
%   For vectors, sumsqr(X) is the sum of the squares of the elements of
%   the vector; i.e., SUM(X.^2).
%
%   For matrices, sumsqr(X) is a row vector containing the application
%   of sumsqr to each column. For N-D arrays, the sumsqr operation
%   is applied to the first non-singleton dimension of X.
%
%   sumsqr(X,DIM) takes the sum along the dimension DIM of X.
%
%   Disciplined convex programming information:
%       If X is real, then sumsqr(X,...) is convex and nonmonotonic in
%       X. If X is complex, then sumsqr(X,...) is neither convex nor
%       concave. Thus, when used in CVX expressions, X must be affine. DIM
%       must be constant.

narginchk(1,2);
y = x .* x;
if nargin == 2
    y = sum( y, dim );
else
    y = sum( y );
end

% Copyright 2005-2016 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.