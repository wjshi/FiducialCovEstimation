function[ ] = PlotInfEig( folderName, InfeigFile, penName )
%Trace plots for leading eigenvalue, ratio of two leading
%eigenvalues, angle between two leading eigenvectors, and condition number
%for covariance matrix estimators.


%addpath ~/Documents/research/matlab/exportfig/
addpath exportfig/    %cluster
load( InfeigFile )
Mtype = {'Eig1' 'Eig1/Eig2' 'Eigvec angle' 'Cond' };

figure('units','normalized','outerposition',[0 0 1 1]);
for type = 1:4
    if type == 1
        M = Eigval1_6;
    elseif type == 2
        M = Eig1Eig2_6;
    elseif type == 3
        M = Theta_6;
    elseif type == 4
        M = Cond_6;
    end
    mtype = Mtype{ type };
    subplot( 4, 1, type  )
    plot( 1:Window, M );
    xlabel( 'Step' )
    title( mtype )              
end

hl = legend( namestr, 'Location', 'south', 'Orientation',...
    'horizontal' );
pos = get(hl,'position');
% Cluster
% pos(2) = pos(2) - 0.2 ; %move the legend down
% pos(1) = pos(1) + 0.1 ; %move the legend to the right
% Local
pos(2) = pos(2) - 0.05 ; %move the legend down

set( hl, 'position', pos ); 
    
% set(gcf,'color',[1 1 1]);
export_fig( sprintf('%s/Infeig_%s.jpg', folderName, penName ) )

% close gcf; clear gcf    

end









    
