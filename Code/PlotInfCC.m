function[ ] = PlotInfCC( folderName, InfFile, InfeigFile, penName )
%Plot confidence curves
addpath exportfig/    %cluster

load( InfFile )
load( InfeigFile )
oracle = length(namestr);

% colors = { 'blue' 'green' 'red' 'cyan' 'magenta' 'yellow' };
colors = { 'blue' 'red' 'yellow' 'magenta' 'green' 'cyan'};
Mtype = { 'GFD' 'Dim' 'D2Sn' 'D2Sig' 'LogD' 'Eig1' 'Eig1/Eig2'...
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
    for j = 1:5
        PlotConfcurve( M( j, : ),  colors{ j } );
        hold on
    end
    plot([ M( oracle ,1 ) M( oracle ,1 ) ],[ 0 1 ] , 'color', colors{ 6 } ) 
    title( Mtype{ ind } );
    hold off
end 

hl = legend( namestr, 'Location', 'south', 'Orientation',...
    'horizontal' );
pos = get(hl,'position');
%local
pos(2) = pos(2) - 0.05 ; %move the legend down
pos(1) = pos(1) - 0.28 ; %move the legend to the left

%cluster
% pos(2) = pos(2) - 0.2 ; %move the legend down
% pos(1) = pos(1) - 0.1 ; %move the legend to the left

set( hl, 'position', pos ); 
    
% set(gcf,'Color','w');
export_fig(sprintf('%s/InfCC_%s.jpg', folderName, penName ))
 
%close gcf; clear gcf    



















end