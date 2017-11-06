function[ ] = PlotInfKD( folderName, InfFile, InfeigFile, penName )

%Plot kernel density for GFD, dim of A, distance to Sn, distance to Sigma,
%log determinant, leading eigenvalue, ratio of two leading
%eigenvalues, angle between two leading eigenvectors, and condition number
%for covariance matrix estimators.


addpath exportfig/    %cluster

load( InfFile )
load( InfeigFile )
oracle = length(namestr);
% colors = { 'blue' 'green' 'red' 'cyan' 'magenta' 'yellow' };
colors = { 'blue' 'red' 'yellow' 'magenta' 'green' 'cyan'};
Mtype = { 'GFD' 'Dim (shifted)' 'D2Sn' 'D2Sig' 'LogD' 'Eig1' 'Eig1/Eig2'...
    'Eigvec angle' 'Cond' };

figure('units','normalized','outerposition',[0 0 1 1]);
for ind = 1:9
    if ind == 1
        M = GFDs;
    elseif ind == 2
        M = DimofAs;
    elseif ind == 3
        M = DtoSns;
    elseif ind == 4
        M = DtoSigmas;
    elseif ind == 5
        M = Ds;
    elseif ind == 6
        M = Eigval1_6;
    elseif ind == 7
        M = Eig1Eig2_6;
    elseif ind == 8
        M = Theta_6;
    else
        M = Cond_6;
    end
    
    subplot( 3, 3, ind )    
    if ind == 2        
        dimAs = M( 1:5, : );
        nbin = max( dimAs(:) ) - min( dimAs(:) );
        [ f, x ] = hist( dimAs', nbin );
        for j = 1:5
            xj = x - 0.6 + 0.2 * j; 
            bar( xj, f( :, j ), 0.15, colors{ j },...
                'LineStyle', 'none');
            hold on
        end
        y6 = 1.1 * max( f(:) );
    else 
        f_max = zeros( 5, 1 );
        for j = 1:5
            [ f, x ] = ksdensity( M( j, : ) );
            f_max( j ) = max( f );
            plot( x, f, 'color', colors{ j } )
            hold on
        end
        y6 = 1.1 * max( f_max );
    end
    
    plot([ M( oracle ,1 ) M( oracle ,1 ) ],[ 0 y6 ] , 'color', colors{ 6 } ) 
    title( Mtype{ ind } );
    hold off
end 

hl = legend( namestr, 'Location', 'south', 'Orientation',...
    'horizontal' );
pos = get(hl,'position');
pos(2) = pos(2) - 0.05 ; %move the legend down
pos(1) = pos(1) - 0.28 ; %move the legend to the left
set( hl, 'position', pos );     
% set(gcf,'color',[1 1 1]);

export_fig( sprintf('%s/InfKD_%s.jpg', folderName, penName ) ) 
% close gcf; clear gcf    



















end