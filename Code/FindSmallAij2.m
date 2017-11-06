function [ optim, optimval ] = FindSmallAij2( Sn, A, n, ith, jth )
%Optimize the normal likelihood when varying the ijth entry of A.

    [ optim, optimval ] = fminsearch( @( x ) Snorm( x, Sn, A, n, ith,...
        jth ), 0 );
    optimval = -optimval;
end

function value = Snorm( x, Sn, A, n, ith, jth )
    Atemp = A;
    Atemp( ith, jth ) = x;
    Stemp = Atemp * Atemp';

    value = n * log( abs( det( Atemp ) ) ) + trace( n * Sn / Stemp ) / 2;
 
end

