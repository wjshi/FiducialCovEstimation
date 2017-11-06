function [ DtoSns, DtoSigmas, DimofAs, Ds, Ts, Es, GFDs, namestr ] = ...
    Inference( dataFile, folderName, mcmcFiles, saveFile, Window, Pen ) 


    load( dataFile ) %e.g. SimData/p4n10/dup1.mat


    nonzeros0 = getnon0s( A );
    logjacob0 = LogJacob( V, A, nonzeros0, 100 );
    FD = logjacob0 + logNormal( V, A );
    if ( nargin == 6 )
        gfd = FD - feval( Pen, A, n ); 
    else 
        gfd = FD;
    end
 
    DimofA0 = sum( sum (A ~= 0 ) );
    DtoSn0 = ForstnerMetric( Sigma, Sn );

    p = size( A, 1 );
    files = dir( mcmcFiles );
    filenum = size( files, 1 );
    GFDs = ones( filenum + 1, Window ) * gfd;
    
    
    
    
    DtoSns = ones( filenum + 1, Window ) * DtoSn0;
    DtoSigmas = DtoSns;
    DimofAs = ones( filenum + 1, Window) * DimofA0;
    Ds = ones( filenum + 1, Window ) * log( abs( det( A ) ) ) * 2;
    Ts = ones( filenum + 1, Window ) * trace( Sigma );
    Es = ones( filenum + 1, Window ) * norm( Sigma );
    namestr = cell( 1, filenum );
    namestr( filenum + 1 ) = cellstr( 'oracle/Sn' );

    
    for i = 1 : filenum
        eval( [ 'load ' folderName '/' files( i ).name ] )
        namestr( i ) = cellstr( files( i ).name( 5:8 ) );
        As = full( SparseAs );

        for scan = 1 : Window
            simA = As( :, ( 1 : p ) + ( scan - 1 ) * p );
            simS = simA * simA';
            DtoSns( i, scan ) = ForstnerMetric( Sn, simS );
            DtoSigmas( i, scan ) = ForstnerMetric( Sigma, simS );
            DimofAs ( i, scan ) = sum( sum( simA ~= 0 ) );
            Ds( i, scan ) = log( abs( det( simA ) ) ) * 2;
            Ts( i, scan ) = trace( simS );
            Es( i, scan ) = norm( simS );
        end
        
        d1 = size( GFD, 1 );
        if d1 == 1
            GFDs( i, : ) = GFD;
        else 
            GFDs( i, : ) = GFD';
        end
        
    end


save ( saveFile, 'DtoSns', 'DtoSigmas', 'DimofAs', 'Ds', 'Ts', 'Es',...
    'GFDs', 'namestr', 'Window' ) 
end


% 
% 
% 
% 
% 
% 
