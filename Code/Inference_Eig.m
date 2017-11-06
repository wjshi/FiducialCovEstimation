function [ Eigval1, Eigval2, Eigvec1, Condnum, Eigval1_6,...
    Eig1Eig2_6, Theta_6, Cond_6 ] = ...
    Inference_Eig( dataFile, folderName, mcmcFiles, saveFile, Window ) 
%%%%%%%%%%%%%%%%%%%%%%%%%%
%Eigval1: the largest eigenvalue
%Eigval2: the second largest eigenvalue
%Eigvec1: the leading eigenvector
%Condum: condition number for the estimated covariance matrices
%Eigval1_6: the largest eigenvalue along with oracle
%Eig1/Eig2_6: the ratio of the largest two eigenvalues along with oracle
%Theta_6: the angle between the leading two eigenvectors along with oracle
%Cond_6: condition indices along with oracle
%%%%%%%%%%%%%%%%%%%%%%%%%%

    load( dataFile ) %e.g. SimData/p4n10/dup1.mat

    files = dir( mcmcFiles );
    filenum = size( files, 1 );    
    
    Eigval1 = ones( filenum, Window );
    Eigval2 = Eigval1;
    Eigvec1 = ones( filenum, p, Window );
    Condnum = ones( filenum, Window );
           
    namestr = cell( 1, filenum );
    namestr( filenum + 1 ) = cellstr( 'oracle/Sn' );
    
    for i = 1 : filenum
        eval( [ 'load ' folderName '/' files( i ).name ] )
        namestr( i ) = cellstr( files( i ).name( 5:8 ) );
        As = full( SparseAs );

        for scan = 1 : Window
            simA = As( :, ( 1 : p ) + ( scan - 1 ) * p );
            simS = simA * simA';
            [ Ve, D ] = eigs( simS, 2 );
            Eigval1( i, scan ) = D( 1, 1 );
            Eigval2( i, scan ) = D( 2, 2 );
            Eigvec1( i, :, scan ) = Ve( :, 1 );
            Condnum( i, scan ) = cond( simS );
        end
        
    end

    
    [ Ve, D ] = eigs( Sigma, 2 );
    eigval1_S = D( 1, 1 );
    eigval2_S = D( 2, 2 );
    eigvec1_S = Ve( :, 1 );
    condnum_S = cond( Sigma ); 

    %Mtype = {'Eig1' 'Eig1/Eig2' 'Eigvec angle' 'Cond' };
    

    M = ones( filenum + 1, Window );
    
    Eigval1_6 = M * eigval1_S;
    Eigval1_6( 1 : filenum, : ) = Eigval1;
    
    Eig1Eig2_6 = M * eigval1_S / eigval2_S; 
    for chain = 1 : filenum
        for scan = 1 : Window
            Eig1Eig2_6( chain, scan ) = Eigval1( chain, scan ) /...
                Eigval2( chain, scan );
        end 
    end
    
    [ a, ~ ] = eigs( Sn, 1 );
    b = eigvec1_S;
    Theta_6 = M * acosd( abs( dot( a, b ) ) );
    for chain = 1 : filenum
        for scan = 1 : Window
            a = Eigvec1( chain, :, scan );
            Theta_6( chain, scan ) = acosd( abs( dot( a, b ) ) );
        end 
    end
    Theta_6 = real(Theta_6); 
    
    Cond_6 = M * condnum_S;
    Cond_6( 1 : filenum, : ) = Condnum;
    

  


save( saveFile, 'Eigval1', 'Eigval2', 'Eigvec1', 'Condnum', 'namestr',...
     'Eigval1_6', 'Eig1Eig2_6', 'Theta_6', 'Cond_6', 'Window') 
end

