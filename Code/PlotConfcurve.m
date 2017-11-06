function[ f, x ] = PlotConfcurve( X, Color )

[ f, x ] = ecdf( X );

for ind = 1: length( f )
    if f( ind ) >= 0.5  
        f( ind ) = 1 - f( ind );
    end
end
f = f*2;


plot( x, f, 'color', Color )

end