function[penalty] = sumPi2( Covariate, n )
%Penaly with form sum(Pi^2/2*log(n))
%P: a vector of nonzeros counts in each row of the covariate matrix. Pi is
%the number of nonzeros in the ith row of covariate matrix.
%n: number of observations.
    P = sum( Covariate' ~= 0 );
    p = length( P );

    penalty = sum( P.^2 / 2 * log( n * p ));
end
