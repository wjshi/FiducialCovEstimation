function [d_AB]=ForstnerMetric(A,B)

%A metric for covariance matrices developed by Wolfgang Forstner and
%Boudewijn Moonen.
%Input:
%   A,B: symmetric, positive definite matrices.
%Output:
%   d_AB: distance metric. d_AB=sqrt(sum(ln(eigenvalues of A*inv(B)).^2))
%Comment:
%   1. d_AB is invariant with respect to affine transformations of the
%   coordinate system.
%   2. d_AB is invariant with respect to an inversion of the matrices.


    d_AB=sqrt(sum(log(eig(A/B)).^2));
end