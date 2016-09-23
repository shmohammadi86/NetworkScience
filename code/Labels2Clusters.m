function [Clusters]=Labels2Clusters(labels)
    UC = unique(labels);
    
    UC = setdiff(UC, -1); % use '-1' as unasigned label 
    
    n = numel(UC);
    Clusters = cell(n, 1);
    for c = 1:n
        c_id = UC(c);
        Clusters{c} = find(labels == c_id);
    end
end
