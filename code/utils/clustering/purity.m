function [ p ] = purity( lab1, lab2 )
    N = numel(lab1);
    
    p = 0;
    UC = unique(lab1);
    for i = 1:numel(UC)
        idx = lab1 == UC(i);
        
        M = 0;
        labels = lab2(idx);
        UC2 = unique(labels);
        for j = 1:numel(UC2)
            count = nnz(labels == UC2(j));
            if(M < count)
                M = count;
            end
        end
        
        p = p + M;
    end
        
    p = p / N;
end

