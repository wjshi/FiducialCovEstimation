function[ f, x ] = PlotKDorHist( X, Color, Isctns, Binnum )

if Isctns == 1
    [ f, x ] = ksdensity( X );
else
    [ f, x ] = hist( X, Binnum );
end



plot( x, f, 'color', Color )

end