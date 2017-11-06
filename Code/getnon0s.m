function[NonZeros] = getnon0s(M)
%Extract rowwise nonzero list for a given matrix M

D=size(M);
r=D(1);
NonZeros=cell(1,r);
for i=1:r
    NonZeros{i}=(find(M(i,:)~=0))';
end

    
