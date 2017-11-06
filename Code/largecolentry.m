function[ A ] = largecolentry( S, MaxC )
%
%Input (S,MaxC)
%   S: Sample covariance matrix.
%   MaxC: Maximum number of non-fixed zeros allowed per column.
%
%Output [A]
%   A: Proposed square-root approximation of S with column constriant.

    A = S^.5;
    C = size( S, 2 );
    q = ( C - MaxC ) / C;
    Q = quantile( A, q, 2 );
    for c = 1 : C
        non0 = gt( A( :, c ), Q( c ) );
        if isequal( non0( c ), 0 )
            disp( 'Diagonal entry of Sn is too small' )
            non0 = zeros( C, 1 );
            non0( c ) = 1;
        end
        A(:,c)=A(:,c).*non0;
    end
        