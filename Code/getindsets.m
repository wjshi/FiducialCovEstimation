function[IndSets] = getindsets(M, MaxC)
%Extract lists of update, birth, death set for a given matrix M
%MaxC is the maximum nonzeros allowed per column


[r,c]=size(M);
IndSets=cell(1,3);
Diag=[1:r; 1:r]';
for j=1:c
    non0sj=find(M(:,j)~=0);
    Lj=length(non0sj);
    S1=([non0sj'; ones(1,Lj)*j])';
    if Lj==MaxC
        S2=[];
    else
        S2=[setdiff(1:r,non0sj); ones(1,r-Lj)*j]';
    end
    S3=setdiff(S1,Diag,'rows');
    IndSets{1}=vertcat(IndSets{1},S1);
    IndSets{2}=vertcat(IndSets{2},S2);
    IndSets{3}=vertcat(IndSets{3},S3);
end