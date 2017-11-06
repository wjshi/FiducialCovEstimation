function [ SparseAs, GFD, a0, b0 ]...
    = mcmcf(initialA, sampleFile, saveFile, burnin, Window,...
    IndexSets, Nonzeros)


load( sampleFile ); %sim11.mat

A0 = initialA;

numblocksample = 100;
s1 = s;
[logjacob, blocksample] = LogJacob(V, A0, Nonzeros, numblocksample);

a0 = 0;
b0 = 0;

w = 1;
As = zeros( p, Window * p );
GFD = ones( 1, Window );
W = burnin + Window;

while w <= W
    [ A0, a0, b0, Nonzeros, blocksample, IndexSets, logjacob ] =...
        Sp1step_0( V, A0, s1, a0, b0, IndexSets, logjacob, Nonzeros,...
        numblocksample, blocksample);
    
    if w > burnin
        wb = w - burnin;
        As( :, ( 1 : p ) + p * (wb - 1) ) = sparse(A0);
        GFD( wb ) = logjacob + logNormal( V, A0 );
    end
    if mod( w, 5000 ) == 0
        disp( w );
%          SparseAs = sparse( As ); %As = full( SparseAs )
%         saveFileNoMat = saveFile(1:end-4);
%         tempFile = sprintf('%s-iter%d.mat', saveFileNoMat, w );
%         save( tempFile, 'initialA', 'SparseAs', 'GFD', 'a0', 'b0', '-v7.3');
    save( saveFile, 'initialA', 'As', 'GFD', 'a0', 'b0',...
        'w', 'burnin', 'Window', '-v7.3'); 
    end
    w=w+1;
end
 SparseAs = As;
 save( saveFile, 'initialA', 'SparseAs', 'GFD', 'a0', 'b0',...
        'burnin', 'Window', '-v7.3');




