function [ penalty ] = MDL( Covariate, n )
%Penaly with form: 
%sum( Pi^2 / 2 * log( n*p ) ) + sum(  log( nchoosek(p, Pi ) ) )
%
%P: a vector of nonzeros counts in each row of the covariate matrix. Pi is
%the number of nonzeros in the ith row of covariate matrix.
%n: number of observations.

    P = sum( Covariate' ~= 0 );
    p = length( P );
    penalty = 0;

    for i = 1 : p
        pi = P( i );
        peni = pi / 2 * log( n * p )+ gammaln( p + 1 ) - ...
            gammaln( p - pi +1 ) - gammaln( pi + 1 ) ;
        penalty = penalty + peni; 
    end





end