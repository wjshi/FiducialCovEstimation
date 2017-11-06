function [  ] = PlotInf( folderName, InfFile, penName )

%Trace plots of GFD, dim of A, distance to Sn, distance to Sigma,
%log determinant, leading eigenvalue.

 addpath exportfig/ 
% addpath ~/Documents/research/matlab/exportfig/
load( InfFile ) %Inf.mat
Mtypes = { 'GFD' 'Dim' 'D2Sn' 'D2Sig' 'LogD' 'Eig1'};
    

figure('units','normalized','outerposition',[0 0 1 1]);
    for Mtype = 1:6
        subplot( 3, 2, Mtype )   
        if Mtype == 1
            M = GFDs;
        elseif Mtype == 2
            M = DimofAs;
        elseif Mtype == 3
            M = DtoSns;
        elseif Mtype == 4
            M = DtoSigmas;
        elseif Mtype == 5
            M = Ds;
        elseif Mtype == 6
            M = Es;
        else
            disp( 'Unknown Mtype.' ) 
        end

        plot( 1:Window, M );
        xlabel ( 'Step' )
        title ( Mtypes(Mtype) )
        hold on          
    end

    hl = legend( namestr, 'Location',...
        'South', 'Orientation', 'horizontal' );
    pos = get(hl,'position');
    pos(1) = pos(1)-0.22; %move the legend to the left by 0.05.
    pos(2) = pos(2)-0.05; %move the legend down by 0.13.
    set(hl,'position',pos);
  hold off

%      set(gcf,'color',[1 1 1])
     export_fig( sprintf('%s/Inf_%s.jpg', folderName, penName ) ) 
%      close gcf; clear gcf  

end



    
