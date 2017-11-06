
function [ A1, a1,  b1, newNonzeros, newBlocksample, newIndexSets,...
    newlogJacob]...
= Sp1step_0( V, A0, s1, a0, b0, IndexSets, logJacob, Nonzeros,...
    numblocksample, Blocksample)
%Approximate one step Metroplis sampling with sparsity restriction on A1 
%when dimension of A0 and sample size are large. 
%%%%%%%%%%%%%%%%%%%%%%%%
%Strucure of A is known%
%%%%%%%%%%%%%%%%%%%%%%%%
%
%Input: ( V, A0, s1, a0, b0, IndexSets, logJacob, Nonzeros, N,...
% numblocksample, Blocksample, Jfcn )
%   V: stacked observations. V=[Y1,...,Yn]'.
%   A0: p*p coefficient matrix A0. A0 needs to have full rank.
%   s1: relatively large entry variance in the "update" move.
%   a0: 3*1 vector. Record numbers of proposed update acceptance in each 
%       move for this step. 
%   b0: 3*1 vector. Record numbers of proposed moves before this step. 
%   IndexSets: a three-part-list consists of possible update indecies, 
%        possible birth indecies and possible death indecies. 
%   logJacob: log-transformed adjusted Jacobian of A0
%   Nonzeros: Nonzeros entries of A0 recorded is a list of column numbers
%       for each row. 
%   numblocksample: number of block determinant calculated. It is only used
%       when dimension of A0 and Y are large, i.e. when exact Jacobian is
%       expensive to compute. 
%   Blocksample: specified block samples for the Jacobian calculation.
%
%Output: [ A1, a1,  b1, newNonzeros, newBlocksample, newIndexSets,...
%     newlogJacob]
%   A1: updated coefficient matrix.
%   a1: Record the number of proposed update acceptance in each 
%       move after this step. In one step case each tuple is either 0 or 1.
%   b1: Record the number of proposed moves after this step. In 
%       one step case each tuple is either 0 or 1. 
%   newNonzeros: Nonzeros entries of A1 recorded is a list of column numbers
%       for each row. 
%   newBlocksample: block samples used for the Jacobian calculation.
%   newIndexSets: a three-part-list consists of possible update indecies, 
%        possible birth indecies and possible death indecies after this step. 
%   newlogJacob: log-transformed adjusted Jacobian of A1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n = size( V, 1 );
    % p = size( A0, 1 );
    Sn = V' * V / n;
    sig0 = A0 * A0';
    a1 = a0;
    b1 = b0;
    newIndexSets = IndexSets;
    newNonzeros = Nonzeros;
    lognorm0 = -n * log( abs( det( A0 ) ) ) - trace( Sn / sig0 ) * n / 2;
    rcond_cutoff = .000001;



%compute acceptance ratio a
    b1( 1 ) = b0( 1 ) + 1;
    %Choose an entry from update index set (UIS)
    L = IndexSets{ 1 };
    ind = L( randi( size( L, 1 ) ), : );
    row = ind( 1 );
    col = ind( 2 );
    Aij = A0( row, col );
    %Sample the new entry from Norm(Aij,s1)
    new_Aij = normrnd( Aij, s1 );
    A1 = A0;
    A1( row, col ) = new_Aij;
    while rcond( A1 ) < rcond_cutoff
        new_Aij = normrnd( Aij, s1 ); 
        A1( row, col ) = new_Aij;
    end
    %Compute log(acceptance ratio)
    [ newlogJacob, newBlocksample ] = LogJacob(V, A1, Nonzeros,...
        numblocksample, Blocksample);
    logJ_rate = newlogJacob - logJacob;
    lognorm1 = -n * log( abs( det( A1 ) ) ) - trace( Sn / ( A1 *...
        A1' ) ) * n / 2;
    log_f_rate = lognorm1 - lognorm0;
    loga = log_f_rate + logJ_rate;
    u2 = rand( 1 );
    if log( u2 ) < loga %Accept A1
         a1( 1 ) = a0( 1 ) + 1;
    else
        %Reject the proposal
        A1 = A0;
        newBlocksample = Blocksample;
        newlogJacob = logJacob;
    end 


end




