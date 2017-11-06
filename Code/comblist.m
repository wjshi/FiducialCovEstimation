function [combs] = comblist(n,Maxpi)
%Create combinations of sample indices for sample size 1, 2, ..., Maxpi.
%n: number of observations.
%Maxpi: maximum sample size.
combs=cell(Maxpi,1);
for i=1:Maxpi
    combs{i}=combntns((1:n),i);
end