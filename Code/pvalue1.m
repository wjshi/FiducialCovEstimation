function [ p1 ] = pvalue1( A0, As, x )

%   input: ( A0, As, x )
%   A0: a p-by-p true covariate matrix
%   As: a p-by-Np matrix consists of N estimated covariate matrices
%   x: an N-by-1 vector indicating direction
%
%Output: [ p1 ]
%   p1: p-value


r0 = x' * ( A0 * A0' ) * x;
p = size( A0, 1 );
N = length( As ) / p;
R_Y = ones( N, 1 );
for scan = 1:N
    simA = As( :, ( 1 : p ) + ( scan - 1 ) * p );
    simS = simA * simA';
    rscan = x' * simS * x;
    if rscan > r0
        R_Y( scan ) = 0;
    end 
end

p1 = mean( R_Y );


end