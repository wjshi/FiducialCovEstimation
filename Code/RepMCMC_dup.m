function[ ] = RepMCMC_dup( p, n, duplicate, burnin, Window )

%Run rjmcmcf for simulated data
%Two penalties considered: MDL & EBIC

configCluster
%configCluster( 'killdevil' )


folder_pn = sprintf( 'SimData/p%dn%d', p, n );
folder_dup = sprintf('%s/dup%dburn%dk', folder_pn, duplicate, burnin/1000);
dataFile = sprintf( 'SimData/Simp%dn%d.mat', p, n );
mkdir( folder_dup ); 
load( dataFile );


rng(n + p*10 + duplicate);
longZ = normrnd( 0, 1, n * p, 1 );
B = kron( eye( n ), A );
Y = B * longZ;
V = zeros( n, p );
for ind = 1 : n
    V( ind, : ) = Y( ( ( ind - 1 ) * p + 1 ) : ( ind * p ) );
end
Sn = V'*V/n;
dupFile = sprintf('%s/dup%d.mat', folder_dup, duplicate);
save(dupFile,'p','s','MaxC','A','non0s','indsets','V','Sn',...
    'Sigma','n', 'duplicate')





A = A; Sn = Sn; MaxC = MaxC;
cc = p - MaxC + 1;

%matlabpool open killdevil32 5
matlabpool open 4
parfor iCase = 1 : 4

    switch iCase
        case 1
            initialA = diag( diag( chol( Sn, 'lower' ) ) );
            saveFile = sprintf('%s/MCMCdcho_MDL.mat', folder_dup );
            rjmcmcf( initialA, dupFile, saveFile, burnin, Window, @MDL );
        case 2
            initialA = largecolentry( Sn^.5, MaxC );
            saveFile = sprintf('%s/MCMCMaxC_MDL.mat', folder_dup );
            rjmcmcf( initialA, dupFile, saveFile, burnin, Window, @MDL );
        case 3
            initialA = chol( Sn, 'lower' );
            if MaxC < p
                SetB2 = nchoosek( 1 : cc, 2 );
                zeroout = [ SetB2( :, 2 ) + MaxC - 1, SetB2( :, 1 ) ];
                zerooutind = ( zeroout( :, 2 ) - 1 ) * p +...
                    zeroout( :, 1 );
                initialA( zerooutind ) = 0;
            end
            saveFile = sprintf('%s/MCMCchol_MDL.mat', folder_dup );
            rjmcmcf( initialA, dupFile, saveFile, burnin, Window, @MDL );
        case 4
            initialA = diag(diag((Sn)^.5));
            saveFile = sprintf('%s/MCMCdiag_MDL.mat', folder_dup );
            rjmcmcf( initialA, dupFile, saveFile, burnin, Window, @MDL );
%         case 5
%             initialA = A;
%             saveFile = sprintf('%s/MCMCtrue_MDL.mat', folder_dup );
%             rjmcmcf( initialA, dupFile, saveFile, burnin, Window, @MDL );
            
%         case 6
%             initialA = diag( diag( chol( Sn, 'lower' ) ) );
%             saveFile = sprintf('%s/MCMCdcho_EBIC.mat', folder_dup );
%             rjmcmcf( initialA, dupFile, saveFile, burnin, Window, @EBIC );
%         case 7
%             initialA = A;
%             saveFile = sprintf('%s/MCMCtrue_EBIC.mat', folder_dup );
%             rjmcmcf( initialA, dupFile, saveFile, burnin, Window, @EBIC );
%         case 8
%             initialA = chol( Sn, 'lower' );
%             if MaxC < p
%                 SetB2 = nchoosek( 1 : cc, 2 );
%                 zeroout = [ SetB2( :, 2 ) + MaxC - 1, SetB2( :, 1 ) ];
%                 zerooutind = ( zeroout( :, 2 ) - 1 ) * p +...
%                     zeroout( :, 1 );
%                 initialA( zerooutind ) = 0;
%             end
%             saveFile = sprintf('%s/MCMCchol_EBIC.mat', folder_dup );
%             rjmcmcf( initialA, dupFile, saveFile, burnin, Window, @EBIC );
%         case 9
%             initialA = diag(diag((Sn)^.5));
%             saveFile = sprintf('%s/MCMCdiag_EBIC.mat', folder_dup );
%             rjmcmcf( initialA, dupFile, saveFile, burnin, Window, @EBIC );
%         case 10
%             initialA = largecolentry( Sn^.5, MaxC );
%             saveFile = sprintf('%s/MCMCMaxC_EBIC.mat', folder_dup );
%             rjmcmcf( initialA, dupFile, saveFile, burnin, Window, @EBIC );
        otherwise
            disp( 'Unkown initial condition.' );
    end
    
end
matlabpool close




%Inference
%MDL
mcmcFiles_MDL = sprintf( '%s/MCMC*_MDL.mat', folder_dup );
saveFile1 = sprintf( '%s/Inf_f%dk_MDL.mat', folder_dup, Window / 1000 );
Inference( dupFile, folder_dup, mcmcFiles_MDL,...
    saveFile1, Window, @MDL ); 
saveFile2 = sprintf( '%s/Infeig_f%dk_MDL.mat', folder_dup, Window/1000 );
Inference_Eig( dupFile, folder_dup, mcmcFiles_MDL, saveFile2, Window );
PlotInf(folder_dup, saveFile1, 'MDL');
PlotInfEig(folder_dup, saveFile2, 'MDL');
PlotInfCC(folder_dup, saveFile1, saveFile2, 'MDL');
PlotInfKD(folder_dup, saveFile1, saveFile2, 'MDL');


% %EBIC
% mcmcFiles_EBIC = sprintf( '%s/MCMC*_EBIC.mat', folder_dup );
% saveFile1 = sprintf( '%s/Inf_f%dk_EBIC.mat', folder_dup, Window / 1000 );
% Inference( dupFile, folder_dup, mcmcFiles_EBIC,...
%     saveFile1, Window, @EBIC ); 
% saveFile2 = sprintf( '%s/Infeig_f%dk_EBIC.mat', folder_dup, Window/1000 );
% Inference_Eig( dupFile, folder_dup, mcmcFiles_EBIC, saveFile2, Window );
% PlotInf(folder_dup, saveFile1, 'EBIC');
% PlotInfEig(folder_dup, saveFile2, 'EBIC');
% PlotInfCC(folder_dup, saveFile1, saveFile2, 'EBIC');
% PlotInfKD(folder_dup, saveFile1, saveFile2, 'EBIC');

end











