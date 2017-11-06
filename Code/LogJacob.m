function [ logJacob, newblocksample ] = LogJacob( V, A, Nonzeros,...
    numblocksample, blocksample )
%Approximate Jacobian when sample size and dimension of A are big.
%Input: (V,A,numblocksample,N,blocksample)
%   V: [Y1',...,Yn']', where Yi's are observations
%   A: a p by p invertible matrix
%   N: When the number of nonzeros in a row of A is greater than N, the 
%       determinates of the submatrices are approximated using sampling 
%       with numblocksample. 
%   Nonzeros: Nonzeros entries of A0 recorded is a list of column numbers
%       for each row. 
%   numblocksample: number of block determinant calculated. It is only used
%       when dimension of A0 and Y are large, i.e. when exact Jacobian is
%       expensive to compute. 
%   blocksample: specified block samples. Optional. 
%Output: [log_jacob,newblocksample]
%   Jacob: Jacobian of (Y,A). 
%   newblocksample: block samples used in the calculation.



    
    p = size( A, 1 );
    n = size( V, 1 );
    Z = ( A \ V' )';
    adjlogJ = ones( p, 1 );
    newblocksample = cell( p, 1 );
    for i = 1 : p
        nonzeroind = Nonzeros{ i };
        pi = length( nonzeroind );
        Z_part = Z( :, nonzeroind );
        if ( ( nargin == 4 ) || ( isempty( blocksample ) ) )
           if ( pi > 1 || n > 50 ) 
               blocks = numblocksample;
               
               blocksamplei = ones( n , pi );
               for row = 1 : blocks 
                    blocksamplei( row, : ) = randperm( n, pi ); 
               end     
                    
           else
               blocks = n;
               blocksamplei = (1 : n)';
           end
        else 
           blocksamplei = blocksample{ i };
           blocks = length( blocksamplei );
        end
        
        blockidet = ones( blocks, 1 );
        for j = 1 : blocks
            detij = det( Z_part( blocksamplei( j, : ), : ) );
            blockidet( j ) = abs( detij );
        end
        newblocksample{ i } = blocksamplei;
        avglogJi = log( sum( blockidet ) ) - log( blocks ); 
        adjlogJ( i ) = avglogJi - pi / 2 * log( 2 )... %1/2 for E(Gaussian)
            - gammaln( ( pi + 1 ) / 2 ); 
    end
    logJacob = sum( adjlogJ );    
 end

