function [ A, n, V, Sn, non0s, indsets, Sigma, s ] = Datagen( ...
    p, n, MaxC, Mean, Std, SampleFile)

    rcond_cutoff = .000001;
    
    A = randsqM( p, MaxC, Mean, Std );
    Sigma = A*A';
    Sig_eig = eigs(Sigma, 2);
    eig1eig2 = abs(Sig_eig(1)/Sig_eig(2));
    A_update = find(A ~= 0);
    newupdate_ind = randsample(A_update, 1);
    while rcond(A)<rcond_cutoff || eig1eig2<2
%         newupdate_ind = randsample(A_update, 1);
        A(newupdate_ind) = A(newupdate_ind) + Std;
        Sigma = A*A';
        Sig_eig = eigs(Sigma, 2);
        eig1eig2 = abs(Sig_eig(1)/Sig_eig(2)); disp(eig1eig2);
    end
    
    non0s = getnon0s( A );
    indsets = getindsets( A, MaxC );

    longZ = normrnd( 0, 1, n * p, 1 );
    B = kron( eye( n ), A );
    Y = B * longZ;
    V = zeros( n, p );
    for i = 1 : n
        V( i, : ) = Y( ( ( i - 1 ) * p + 1 ) : ( i * p ) );
    end

    Sn = V'*V/n;
    s = Std;

save(SampleFile,'p','s','A','eig1eig2', 'MaxC','non0s','indsets',...
    'V','Sn', 'Sigma','n')