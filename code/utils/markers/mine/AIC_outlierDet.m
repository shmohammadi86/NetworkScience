function [ feature_count ] = AIC_outlierDet( v, rel_prob_threshold)
    % Preference is given to the smaller selection of top-ranked elements
    % if bigger selection(s) are not different enough (in terms of relative probability) 
    if(~isvector(v))
        error('v should be a vector');
    elseif(std(v) == 0)
        feature_count = numel(v);
        return;
    end
    
    if(nargin == 1)
        rel_prob_threshold = 1;
    end
  
    m = min(v); M = max(v);    
    m_count = nnz(v== m); M_count = nnz(v == M);
    if( m_count + M_count == numel(v) )
        feature_count = M_count;
        return;
    end        
        
    z = (v - mean(v)) ./ std(v) ;
    [~, perm] = sort(z, 'descend');    
    
    min_AIC = Inf;
    feature_count = 1;    
    separated = false;
    
    l = numel(z);
     % backward sweep
    for s = l-m_count-1:-1:M_count
        n = l - s;
        logStd = log10(std(z(perm(s+1:end))));
        AIC = 2*(n*logStd + (sqrt(2)*s*logfactorial(n)/n) );
        rel_prob = exp( (AIC - min_AIC) / 2 );
        if( rel_prob < rel_prob_threshold )
            min_AIC = AIC;
            feature_count = s;
        end
    end 
    
    % forward sweep
    for n = M_count+1:l-m_count
        s = l - n;
        logStd = log10(std(z(perm(1:n))));
        AIC = 2*(n*logStd + (sqrt(2)*s*logfactorial(n)/n) );
        rel_prob = exp( (AIC - min_AIC) / 2 );
        if( rel_prob < rel_prob_threshold )
            min_AIC = AIC;
            feature_count = n;
        end
    end    
end

