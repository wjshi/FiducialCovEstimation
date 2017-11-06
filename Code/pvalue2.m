function [ p2 ] = pvalue2( A0, As, V_obs, x )

%   input: ( A0, As, V_obs, x )
%   A0: a p-by-p true covariate matrix
%   As: a p-by-Np matrix, consisting of N estimated covariate matrices
%   V_obs: an N-by-1 vector of observations 
%   x: an N-by-1 vector indicating direction
%
%Output: [ p2 ]
%   p2: p-value


fp0 = flatpvalue( A0, V_obs, x );
p = size( A0, 1 );
N = length( As ) / p;
R_Y = ones( N, 1 );

for scan = 1:N
    simA = As( :, ( 1 : p ) + ( scan - 1 ) * p );
    fpcan = flatpvalue( simA, V_obs, x );
    if fpcan > fp0
        R_Y( scan ) = 0;
    end 
end

p2 = mean( R_Y );


end