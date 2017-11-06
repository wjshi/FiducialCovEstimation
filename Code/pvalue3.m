function [ p3 ] = pvalue3( As, V_obs, x )

%   input: ( As, V_obs, x )
%   A0: a p-by-p true covariate matrix
%   As: a p-by-Np matrix, consisting of N estimated 
%   covariate matrices
%   V_obs: an N-by-1 vector of observations 
%   x: an N-by-1 vector indicating direction
%
%Output: [ p3 ]
%   p3: p-value

p = size( As, 1 );
N = length( As ) / p;
S_Y = ones( N, 1 );

for scan = 1:N
    simA = As( :, ( 1 : p ) + ( scan - 1 ) * p );
    S_Y( scan ) = flatpvalue( simA, V_obs, x );
end

p3 = mean( S_Y );


end