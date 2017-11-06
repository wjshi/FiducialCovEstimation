 function [ logNormf ] = logNormal( V,A )
%Log transformation of normal density for stacked multinormal observations, Y, given coefficient matrix A. 
%
%input: p*p coefficient matrix A, and stacked
%observation matrix V with dimension n*p. 
%A needs to have full rank.
%output: GFD of (Y,A) without normalizing constant.

    p=length(A);
    n=length(V);
    Sn=V'*V/n;
    sig=A*A';
    logf=-n*p/2*log(2*pi)-n*log(abs(det(A)))-trace(n*Sn/sig)/2;
    logNormf=logf;
end
