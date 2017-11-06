function [ A1, a1,  b1, newNonzeros, newBlocksample, newIndexSets,...
    newlogJacob]...
= SpRJ1step_1( V, A0, s1, a0, b0, IndexSets, logJacob, Nonzeros,...
    numblocksample, Blocksample, Pmoves, MaxC, MDL )
%Approximate one step reversible jump MCMC with sparsity restriction on A1 
%when dimension of A0 and sample size are large. 
%Accept birth/death proposals according to likelihood (no Jacobian)
%%%%%%%%%%%%%%%%%%%%%%%
%No constraint assumed%
%%%%%%%%%%%%%%%%%%%%%%%
%
%Three possible moves: update, birth, and death.
%
%Input: ( V, A0, s1, a0, b0, IndexSets, logJacob, Nonzeros, N,...
% numblocksample, Blocksample, Pmoves, MaxC, MDL )
%   V: stacked observations. V=[Y1,...,Yn]'.
%   A0: p*p coefficient matrix A0. A0 needs to have full rank.
%   s1: relatively large entry variance in the "update" move.
%   a0: 3*1 vector. Record numbers of proposed update acceptance in each 
%       move for this step. 
%   b0: 3*1 vector. Record numbers of proposed moves before this step. 
%   IndexSets: a three-part-list consists of possible update indecies, 
%        possible birth indecies and possible death indecies. 
%   logJacob: log-transformed adjusted Jacobian of A0
%   Nonzeros: Nonzeros entries of A0 recorded is a list of column numbers
%       for each row. 
%   numblocksample: number of block determinant calculated. It is only used
%       when dimension of A0 and Y are large, i.e. when exact Jacobian is
%       expensive to compute. 
%   Blocksample: specified block samples for the Jacobian calculation.
%   Pmoves: 3*3 vector. Probabilities of proposing each move update, birth,
%       and death.
%   MDL: Penalty function.
%
%Output: [ A1, a1,  b1, newNonzeros, newBlocksample, newIndexSets,...
%     newlogJacob]
%   A1: updated coefficient matrix.
%   a1: 3*1 vector. Record numbers of proposed update acceptance in each 
%       move after this step. In one step case each tuple is either 0 or 1.
%   b1: 3*1 vector. Record numbers of proposed moves after this step. In 
%       one step case each tuple is either 0 or 1. 
%   newNonzeros: Nonzeros entries of A1 recorded is a list of column numbers
%       for each row. 
%   newBlocksample: block samples used for the Jacobian calculation.
%   newIndexSets: a three-part-list consists of possible update indecies, 
%        possible birth indecies and possible death indecies after this step. 
%   newlogJacob: log-transformed adjusted Jacobian of A1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Choose a move type
    u1 = rand( 1 );
    type = sum( u1 > cumsum( Pmoves( 1, : ) ) ) + 1;

    n = size( V, 1 );
    % p = size( A0, 1 );
    Sn = V' * V / n;
    sig0 = A0 * A0';
    a1 = a0;
    b1 = b0;
    newIndexSets = IndexSets;
    newNonzeros = Nonzeros;
    lognorm0 = -n * log( abs( det( A0 ) ) ) - trace( Sn / sig0 ) * n / 2;

    rcond_cutoff = .000001;



    %compute acceptance ratio a
    if type == 1 %update
        b1( 1 ) = b0( 1 ) + 1;
        %Choose an entry from update index set (UIS)
        L = IndexSets{ 1 };
        ind = L( randi( size( L, 1 ) ), : );
        row = ind( 1 );
        col = ind( 2 );
        Aij = A0( row, col );
        %Sample the new entry from Norm(Aij,s1)
        new_Aij = normrnd( Aij, s1 );
        A1 = A0;
        A1( row, col ) = new_Aij;
        while rcond( A1 ) < rcond_cutoff
            new_Aij = normrnd( Aij, s1 ); 
            A1( row, col ) = new_Aij;
        end
        %Compute log(acceptance ratio)
        [ newlogJacob, newBlocksample ] = LogJacob( V, A1, Nonzeros,...
            numblocksample, Blocksample );
        logJ_rate = newlogJacob - logJacob;
        lognorm1 = -n * log( abs( det( A1 ) ) ) - trace( Sn / ( A1 *...
            A1' ) ) * n / 2;
        log_f_rate = lognorm1 - lognorm0;
        loga = log_f_rate + logJ_rate;
        u2 = rand( 1 );
        if log( u2 ) < loga %Accept A1
             a1( 1 ) = a0( 1 ) + 1;
        else
            %Reject the proposal
            A1 = A0;
            newBlocksample = Blocksample;
            newlogJacob = logJacob;
        end 

    elseif type == 2 %birth
        b1( 2 ) = b0( 2 ) + 1;
        %Choose an entry from birth index set (BIS)
        L = IndexSets{ 2 };
        Lb = size( L, 1 );
        indb = randi( Lb );
        ind = L( indb, : );
        row = ind( 1 );
        col = ind( 2 );
        newNonzeros = Nonzeros;
        newNonzeros{ row } = union( Nonzeros{ row }, col );  
        %Compute optimal Aij that maximizes the normal likelihood
        Aijcenter = FindSmallAij2( Sn, A0, n, row, col );
        A1 = A0;
        A1( row, col ) = Aijcenter;
        %Compute s2 using zeroth order method
        logJacob1 = LogJacob( V, A1, newNonzeros, numblocksample );
        lognorm1 = -n * log( abs( det( A1 ) ) ) - trace( Sn / ( A1 *...
            A1' ) ) * n / 2;
        %Consider reversed move
        L_D = [ row col; IndexSets{ 3 } ];  %Form new DIS
        rows_D = L_D( :, 1 );
        cols_D = L_D( :, 2 );
        %Choose an entry from new DIS
        Ld = size( L_D, 1 );   
        Fdr = ones( 1, Ld ); %Rate ratio of delete an entry in new DIS  
        for k = 1 : Ld
            if k == 1 %Back to A0
    %                 logJdk = logJacob - logJacob1;
                    lognormk = lognorm0 - lognorm1;
            else
                temprow = rows_D( k );
                tempcol = cols_D( k );
                tempA = A1;
                tempA( temprow, tempcol ) = 0;

                if rcond( tempA ) > rcond_cutoff %Non-ill-condition tempA
    %                 tempNon0s_D = newNonzeros;
    %                 tempNon0s_D{ temprow } = setdiff( tempNon0s_D{... 
    %                     temprow }, tempcol );
    %                 logJdk = LogJacob( V, tempA, tempNon0s_D,...
    %                     numblocksample ) - logJacob1;
                    tempSig = tempA * tempA';
                    lognormk = -n * log( abs( det( tempA ) ) ) -...
                        trace( n * Sn / tempSig ) / 2 - lognorm1;
                else 
                    lognormk = -Inf;
                end

            end
    %         Fdr( k ) = exp( logJdk + lognormk );
            Fdr(k) = exp( lognormk );
        end
        fd = Fdr( 1 ) / sum( Fdr ); %Rate of picking added entry to delete
        %Fiducial ratio of A0 over tempA
        tempf_rate = Fdr( 1 ) * exp( logJacob - logJacob1 ); 
        %Move ratio of birth to tempA over death to A0
        tempr_rate = ( 1 / Lb * Pmoves( 1, 2 ) ) / ( fd * Pmoves( 2,...
            3 ) * sqrt( 2 * pi ) ); 
        %Penalty ratio of A0 over tempA
        pen0 = feval( MDL, A0, n);
        pen1 = feval( MDL, A1, n);
        tempPen_rate = exp( -pen0 + pen1);
        %Sampling standard deviation
        s2 = tempf_rate * tempPen_rate * tempr_rate;  
        %Sample new entry from Norm(Aijcenter,s2);
        new_Aij = normrnd( Aijcenter, s2 ); 
        A1( row, col ) = new_Aij; %Proposed A1

        if rcond( A1 ) > rcond_cutoff %Non-ill-condition A1
            [ newlogJacob, newBlocksample ] = LogJacob( V, A1,...
                newNonzeros, numblocksample );
            newlognorm = -n * log( abs( det( A1 ) ) ) - trace( Sn /...
                ( A1 * A1' ) ) * n / 2;
            %Compute likelihood of setting A1(row,col)=0
            for k = 1 : Ld
                if k == 1 
    %                 logJdk = logJacob - newlogJacob;
                    lognormk = lognorm0 - newlognorm;
                else
                    temprow = rows_D( k );
                    tempcol = cols_D( k );
                    tempA = A1;
                    tempA( temprow, tempcol ) = 0;

                    if rcond( tempA ) > rcond_cutoff                    
    %                     tempNon0s_D = newNonzeros;
    %                     tempNon0s_D{ temprow } = setdiff( tempNon0s_D{...
    %                         temprow }, tempcol );
    %                     logJdk = LogJacob( V, tempA, tempNon0s_D,...
    %                         numblocksample ) - newlogJacob;
                        tempSig = tempA * tempA';
                        lognormk = -n * log( abs( det( tempA ) ) ) -...
                            trace( n * Sn / tempSig ) / 2 - newlognorm;
                    else
                        lognormk = -Inf;
                    end

                end
    %             Fdr( k ) = exp( logJdk + lognormk );
                Fdr( k ) = exp( lognormk );
            end
            fd = Fdr( 1 ) / sum( Fdr );
            %Compute log(acceptance ratio)
            logf_rate = -log( Fdr( 1 ) ) + newlogJacob - logJacob;
            logPen_rate = -log( tempPen_rate );
            r_rate = ( fd * Pmoves( 2, 3 ) ) / ( 1 / Lb * Pmoves( 1,...
                2 ) * normpdf( new_Aij, Aijcenter, s2 ) );
            loga = logf_rate + logPen_rate + log( r_rate );

            u2 = rand( 1 );
            if log( u2 ) < loga %Accept the proposal
                a1( 2 ) = a0( 2 ) + 1;     
                %Add (row,col) to UIS and keep new DIS
                newIndexSets{ 1 } = [ row col; IndexSets{ 1 } ]; 
                newIndexSets{ 3 } = L_D;
                %Update BIS
                colC = sum( A1 ( :, col ) ~= 0 );
                if colC == MaxC %At maximum number of nonzeros allowed 
                    %Remove (~,col) from BIS
                    newIndexSets{ 2 } = L( L( :, 2 ) ~= col, : );
                else
                    %Remove (row,col) from BIS
                    newIndexSets{ 2 } = removerows( L, indb );  
                end    
            else %Reject the proposal
                A1 = A0;
                newlogJacob = logJacob;
                newBlocksample = Blocksample;
                newNonzeros = Nonzeros;
            end 

        else
            %Reject the proposal
            A1 = A0;
            newlogJacob = logJacob;
            newBlocksample = Blocksample;
            newNonzeros = Nonzeros;
        end
    else 
        b1( 3 ) = b0( 3 ) + 1;  %death
        %Choose an entry from death index set
        L = IndexSets{ 3 };
        rows_D = L( :, 1 );
        cols_D = L( :, 2 );
        Ld = size( L, 1 );
        Fdr = ones( 1, Ld ); %Compute the likelihood of picking each entry
        for k = 1 : Ld
            temprow = rows_D( k );
            tempcol = cols_D( k );
            tempA = A0;
            tempA( temprow, tempcol ) = 0;

            if rcond( tempA ) > rcond_cutoff

    %             tempNon0s_D = Nonzeros;
    %             tempNon0s_D{ temprow } = setdiff( tempNon0s_D{... 
    %                 temprow }, tempcol );
    %             logJdk = LogJacob( V, tempA, tempNon0s_D,...
    %                 numblocksample ) - logJacob;
                tempSig = tempA * tempA';
                lognormk = -n * log( abs( det( tempA ) ) ) -...
                    trace( n * Sn / tempSig ) / 2 - lognorm0;
    %             Fdr( k ) = exp( logJdk + lognormk );
                Fdr( k ) = exp( lognormk );
            else
                Fdr( k ) = 0;        
            end

        end

        if sum( Fdr ) == 0 %No deletion results in non-ill-condition A1
            %Reject the proposal
            loga = -Inf;
        else        
            fds = Fdr / sum( Fdr );
            u3 = rand( 1 );
            ind = sum( u3 > cumsum( fds ) ) + 1;
            row = rows_D( ind );
            col = cols_D( ind );
            A1 = A0;
            A1( row, col ) = 0; %Set A0(row, col)=0 => porposed A1
            newNonzeros = Nonzeros;
            newNonzeros{ row } = setdiff( Nonzeros{ row }, col );
            [ newlogJacob, newBlocksample ] = LogJacob( V, A1,...
                newNonzeros, numblocksample );
            lognorm1 = -n * log( abs( det( A1 ) ) ) - trace( n * Sn /...
                ( A1 * A1' ) ) / 2;
            indd = find( L( :, 1 ) == row & L( :, 2 ) == col ); 
            %remove (row,col) from death index set
            L_D = removerows( L, indd ); 
            Rows = find( A1( :, col ) ~= 0 );
            if length( Rows ) == MaxC - 1 %A1 is near full in column col
                p = size( A1, 1 );
                %Add all nonzero entries in col to BIS
                addtoL_B = [ setdiff( 1 : p, Rows ); ones( 1,...
                    p - MaxC + 1 ) * col ]';
            else 
                %Add (row, col) to BIS
                addtoL_B = [ row col ];
            end
            L_B = [ addtoL_B; IndexSets{ 2 } ]; 

            %Consider the reverse move A1->A0 (birth)
            Lb = size( L_B, 1 );
            %Compute optimal Aij
            Aijcenter = FindSmallAij2( Sn, A1, n, row, col );
            tempA0 = A1;
            tempA0( row, col ) = Aijcenter;
            %Compute s2
            templogJacob0 = LogJacob( V, tempA0, Nonzeros,...
                numblocksample );
            templognorm0 = -n * log( abs( det( tempA0 ) ) ) - trace(...
                Sn / ( tempA0 * tempA0' ) ) * n / 2;
            rows_D = L( :, 1 );
            cols_D = L( :, 2 );

            %compute the likelihood of choosing an entry from DIS
            for k = 1 : Ld
                if k == 1
        %             logJdk = newlogJacob - templogJacob0;
                    lognormk = lognorm1 - templognorm0;
                else
                    temprow = rows_D( k );
                    tempcol = cols_D( k );
                    tempA = tempA0;
                    tempA( temprow, tempcol ) = 0;
                    tempSig = tempA * tempA';
                    
                    if rcond( tempA ) > rcond_cutoff
        %                 tempNon0s_D = Nonzeros;
        %                 tempNon0s_D{ temprow } = setdiff( tempNon0s_D{...
        %                     temprow }, tempcol );
        %                 logJdk = LogJacob( V, tempA, tempNon0s_D,...
        %                     numblocksample ) - templogJacob0;
                        lognormk = -n * log( abs( det( tempA ) ) ) -...
                            trace( Sn / tempSig ) * n / 2 - templognorm0;
                    else
                        lognormk = -Inf;

                    end
                end
        %             Fdr( k ) = exp( logJdk + lognormk );
                    Fdr( k ) = exp( lognormk );
            end
            fd = Fdr( 1 ) / sum( Fdr );
            tempf_rate = Fdr( 1 ) * exp( newlogJacob - templogJacob0 );
            tempr_rate = ( 1 / Lb * Pmoves( 3, 2 ) ) / ( fd * Pmoves(...
                1, 3 ) * sqrt ( 2 * pi ) ); %compute temp rate of moves
            pen0 = feval( MDL, A0, n);
            pen1 = feval( MDL, A1, n);
            tempPen_rate = exp( -pen1 + pen0 );
            s2 = tempf_rate * tempPen_rate * tempr_rate;  
            %Compute acceptance ratio
            log_f_rate = lognorm1 - lognorm0;
            logJ_rate = newlogJacob - logJacob;
            logPen_rate = log( tempPen_rate );
            r_rate = ( 1 / Lb * Pmoves( 3, 2 ) * normpdf( A0( row, col ),...
                Aijcenter, s2 ) ) / ( fds( ind ) * Pmoves( 1, 3 ) );
            loga = log_f_rate + logJ_rate + logPen_rate + log( r_rate );
        end
        
        u2 = rand( 1 );
        if log( u2 ) < loga %Accept the proposal
            a1( 3 ) = a0( 3 ) + 1;
            Lu = IndexSets{ 1 };    %Remove (row, col) from update set.
            indu = find( Lu( :, 1 ) == row & Lu( :, 2 ) == col );       
            newIndexSets{ 1 } = removerows( Lu, indu );
            newIndexSets{ 2 } = L_B;
            newIndexSets{ 3 } = L_D;
        else     
            %Reject the proposal
            A1 = A0;
            newlogJacob = logJacob;
            newBlocksample = Blocksample;
            newNonzeros = Nonzeros;
        end 

    end








end




