function [ logri, logDi ] = Logri( V, Gi, numblocksample )
%Compute logri
%   V: [Y1; ...; Yn]', n by p
%   Gi: list of tuples in clique Gi
%   numblocksample: number of block determinant calculated. It is only used
%       when nchoosek(n, gi) is large. gi = length(Gi)


    n = size(V,1);
    gi = length(Gi);
    Vi = V(:, Gi);
    logDi = LogD_gi(V, Gi, numblocksample);
    logri = logMvGamma(n/2, gi) + logDi + gi * (gi - 1)/2 *...
        log(2 * pi) - n/2 * log(abs(det(Vi' * Vi / n)));
    
end
    
