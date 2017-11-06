function [ penalty ] = EBIC( Covariate, n )
%Penaly with form: 
%sum( Pi / 2 * log( n ) ) + log( nchoosek( p^2, sum( Pi ) ) ) )
%
%P: a vector of nonzeros counts in each row of the covariate matrix. Pi is
%the number of nonzeros in the ith row of covariate matrix.
%n: number of observations.

    P = sum( Covariate' ~= 0 );
    Dim = sum( P );
    p = length( P );

    ebic = ones( Dim, 1 );
    for s = 1 : Dim
        ebic( s ) = log( ( p^2 - Dim + s ) / s );
    end
    EBIC = sum( ebic );

    penalty = Dim / 2 * log( n * p ) + EBIC;



end