function [ logDi ] = LogD_gi( V, Gi, numblocksample )

%   Compute gi*log(avg(det(V^i_\i)))
%   V: n by p observation matrix
%   Gi: indices for ith clique
%   numblocksample: number of subsamples

    n = size(V, 1);
    Vi = V(:, Gi);
    gi = length(Gi);
        
    if pi > 1 || n > 50  
        sub_det = ones(numblocksample, 1);
        for ind = 1:numblocksample
            subM = Vi(randperm(n, gi), :);
            sub_det(ind) = abs(det(subM));
        end
    else
        sub_ind = nchoosek(1:n, gi);
        sub_num = length(sub_ind);
        sub_det = ones(sub_num, 1);
        for ind = 1:sub_num
            subM = Vi(sub_ind(ind, :), :);
            sub_det(ind) = abs(det(subM));
%             sub_det(ind) = gi * log(abs(det(subM)));
        end
    end
    logDi = gi * log(mean(sub_det));
%     logDi = gi * (log(mean(sub_det)) + gammaln(n + 1) - ...
%         gammaln(gi + 1) - gammaln(n - gi + 1));
    
%     logDi = mean(sub_det);
    
 end
