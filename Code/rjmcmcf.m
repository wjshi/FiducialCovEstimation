function [ SparseAs, GFD, a0, b0 ]...
    = rjmcmcf( initialA, sampleFile, saveFile, burnin, Window, Pen)


load( sampleFile ); %simp4n10.mat

%p = size( initialA, 1 );
%MaxC = 3;
cc = p - MaxC + 1;

A0 = initialA;
nonzeros = getnon0s( A0 );
indexsets = getindsets( A0, MaxC );

numblocksample = 100;
s1 = 10;
[ logjacob, blocksample ] = LogJacob( V, A0, nonzeros, numblocksample );

a0 = zeros( 3, 1 );
b0 = zeros( 3, 1 );

w = 1;
As = zeros( p, Window * p );
GFD = ones( 1, Window );
W = burnin + Window;

while w <= W
    B = indexsets{ 2 };
    Bset = size( B, 1 );
    if Bset( 1 ) == 0 %full
        Pmoves = [ .5 0 .5; 1 1 1; .5 .25 .25 ];
    else
        BB = B( :, 1 );
        Dset = size( indexsets{ 3 }, 1 );
        if Bset(1) == cc && isequal( BB == mean( BB ), cc ) %nearly full
            Pmoves = [ .5 .25 .25; .5 0 .5; .5 .25 .25 ];
        elseif Dset( 1 ) == 0 %diagonal
            Pmoves = [ .5 .5 0; .5 .25 .25; 1 1 1 ];
        elseif Dset( 1 ) == 1 %one death option (nearly diagonal)
            Pmoves = [ .5 .25 .25; .5 .25 .25; .5 .5 0 ];
        else
            Pmoves = [ .5 .25 .25; .5 .25 .25; .5 .25 .25 ];
        end
    end

    [ A0, a0, b0, nonzeros, blocksample, indexsets, logjacob ] =...
        SpRJ1step_1( V, A0, s1, a0, b0, indexsets, logjacob, nonzeros,...
        numblocksample, blocksample, Pmoves, MaxC, Pen );  
 
    
    if w > burnin
        wb = w - burnin;
        As( :, ( 1 : p ) + p * (wb - 1) ) = sparse(A0);
        pen = feval( Pen, A0, n );
        GFD( wb ) = logjacob + logNormal( V, A0 ) - pen;
    end
    if mod( w, 5000 ) == 0
        disp( w );
%          SparseAs = sparse( As ); %As = full( SparseAs )
%         saveFileNoMat = saveFile(1:end-4);
%         tempFile = sprintf('%s-iter%d.mat', saveFileNoMat, w );
%         save( tempFile, 'initialA', 'SparseAs', 'GFD', 'a0', 'b0', '-v7.3');
    save( saveFile, 'initialA', 'As', 'GFD', 'a0', 'b0',...
        'w', 'burnin', '-v7.3');

 
    end
    w=w+1;
 end
SparseAs = As;
 save( saveFile, 'initialA', 'SparseAs', 'GFD', 'a0', 'b0',...
        'burnin', 'Window', '-v7.3');
end



