function[ ] = PlotInfCCsmall( folderName, InfFile, InfeigFile, penName )
%Plot confidence curves
addpath exportfig/    %cluster

load( InfFile )
load( InfeigFile )
oracle = length(namestr);
newnamestr = namestr;
newnamestr(oracle) = cellstr('oracle');
newnamestr(oracle + 1) = cellstr('Sn');

% colors = { 'blue' 'green' 'red' 'cyan' 'magenta' 'yellow' };
colors = { 'blue' 'red' 'yellow' 'magenta' 'green' 'cyan' 'cyan'};
Mtype = { 'GFD' 'D2Sig' 'LogD' 'Eigvec angle' };


% figure('units','normalized','outerposition',[0 0 1 1]);
for ind = 1:4
    if ind == 1
        M = GFDs;
    elseif ind == 2
        M = DtoSigmas;
    elseif ind == 3
        M = Ds;
    elseif ind == 4
        M = Theta_6;
    end
    
    subplot( 2, 2, ind )    
    for j = 1:5
        PlotConfcurve( M( j, : ),  colors{ j } );
        hold on
    end
    if mod(ind, 2) == 0
        plot([M( oracle ,1 ) M( oracle ,1 )],[0 0], '-',...
            'color', colors{ oracle });
        hold on
        plot([ M( oracle ,1 ) M( oracle ,1 ) ],[ 0 1 ] , '-.' , ...
            'color', colors{ oracle + 1 })
    else
        plot([ M( oracle ,1 ) M( oracle ,1 ) ],[ 0 1 ] ,...
            'color', colors{ oracle } )
        hold on 
        plot([M( oracle ,1 ) M( oracle ,1 )],[0 0],'-.',...
            'color', colors{ oracle + 1 });
    end
    title( Mtype{ ind } );
    hold off
end 


hl = legend( newnamestr, 'Location', 'south', 'Orientation',...
    'horizontal' );
pos = get(hl,'position');
%local
pos(2) = pos(2) - 0.1 ; %move the legend down
pos(1) = pos(1) - 0.22 ; %move the legend to the left

%cluster
% pos(2) = pos(2) - 0.2 ; %move the legend down
% pos(1) = pos(1) - 0.1 ; %move the legend to the left

set( hl, 'position', pos ); 
    
% set(gcf,'Color','w');
% export_fig(sprintf('%s/InfCCsmall_%s.jpg', folderName, penName ))
 
%close gcf; clear gcf    



















end