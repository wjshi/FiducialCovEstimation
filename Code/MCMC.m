function[ ] = MCMC( p, n, burnin, Window )

%Run mcmcf for simulated data

%configCluster
%configCluster( 'killdevil' )


folder_pn = sprintf( 'SimData/p%dn%d', p, n );
folderName = sprintf( '%s/MCMCburn%dk', folder_pn, burnin/1000 );
dataFile = sprintf( 'SimData/Simp%dn%d.mat', p, n );
mkdir( folderName ); 
load( dataFile );


mu = 0; s = 10;
A = A; Sn = Sn; non0s = non0s;

indexsets = getindsets( A, MaxC );


%matlabpool open killdevil32 5
% matlabpool open 10
parfor iCase = 1:5
    switch iCase
        case 1
            initialA = diag( diag( chol( Sn, 'lower' ) ) );
            saveFile = sprintf('%s/MCMCdcho_MC.mat', folderName );
            mcmcf(initialA, dataFile, saveFile, burnin, Window,...
                indexsets, non0s);
        case 2
            initialA = A;
            saveFile = sprintf('%s/MCMCtrue_MC.mat', folderName );
            mcmcf(initialA, dataFile, saveFile, burnin, Window,...
                indexsets, non0s);
        case 3
            initialA = (Sn)^.5;
            Zeros = A == 0;
            initialA(Zeros) = 0;
            saveFile = sprintf('%s/MCMCSnPa_MC.mat', folderName );
            mcmcf(initialA, dataFile, saveFile, burnin, Window,...
                indexsets, non0s);
        case 4
            initialA = diag(diag((Sn)^.5));
            saveFile = sprintf('%s/MCMCdiag_MC.mat', folderName );
            mcmcf(initialA, dataFile, saveFile, burnin, Window,...
                indexsets, non0s);
        case 5
            rng(n + p*3 + 7)
            initialA = A;
            NonZeros = A ~= 0;
            initialA(NonZeros) = normrnd(mu, s, sum(sum(NonZeros)), 1);
            saveFile = sprintf('%s/MCMCrand_MC.mat', folderName );
            mcmcf(initialA, dataFile, saveFile, burnin, Window,...
                indexsets, non0s);
         otherwise
            disp( 'Unkown initial condition.' );
    end
    
end
% matlabpool close




%Inference
mcmcFiles_MC = sprintf( '%s/MCMC*_MC.mat', folderName );
saveFile1 = sprintf( '%s/Inf_f%dk_MC.mat', folderName, Window/1000 );
Inference(dataFile, folderName, mcmcFiles_MC,...
    saveFile1, Window); 
saveFile2 = sprintf( '%s/Infeig_f%dk_MC.mat', folderName, Window/1000 );
Inference_Eig( dataFile, folderName, mcmcFiles_MC, saveFile2, Window );

PlotInf(folderName, saveFile1, 'MCMC');
PlotInfEig(folderName, saveFile2, 'MCMC');
PlotInfCC(folderName, saveFile1, saveFile2, 'MCMC');
PlotInfKD(folderName, saveFile1, saveFile2, 'MCMC');

end











