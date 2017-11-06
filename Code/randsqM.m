function [ M ] = randsqM( Ncol, Maxc, Mean, Std)
%Generate random matrix M with column constraint. 
%Input( Ncol, Maxc, Mean, Std)
%   Ncol: number of columns for M
%   Maxc: max number of nonzeros allowed per column
%   Mean & Std: mean & standard deviation of random normal generation of
%   entries for M. 


    M = diag( normrnd( Ncol, Std, Ncol, 1));
%     M = diag( Mean + (Std - Mean) * rand( Ncol, 1 ) );
    
    Nrow = Ncol; 

    for col = 1 : Ncol
        rows = randsample( Maxc - 1, 1 );
        turnon = randsample( Nrow, rows );
        M( turnon, col ) = normrnd( Mean, Std, rows, 1);
%         M( turnon, col ) = Mean + (Std - Mean) * rand( rows, 1 );
    end
    
end