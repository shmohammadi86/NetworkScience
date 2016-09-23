function [ SNN_labs] = run_SNNCliq( expression, varargin )
    params = inputParser;    
    params.addParamValue('k', 3      ,@(x) isscalar(x));  % The size of the nearest neighbor list
    params.addParamValue('distance', 'euclidean'      ,@(x) ischar(x));     
    params.addParamValue('r', 0.7      ,@(x) isscalar(x)); % specifies the density threshold of quasi-cliques. Changing r can fine-tune the granularity of cluster output: increasing r results in granulated clusters; reducing r leads to bigger clusters.
    params.addParamValue('m', 0.5      ,@(x) isscalar(x)); % specifies the threshold on the overlapping rate for merging. Change m has the same effect as changing r.
    
    params.parse(varargin{:});
    par = params.Results;    
    
    if(20 < max(nonzeros(expression))) % if "not log-normalized," ...
        expression = log2(expression + 1);
    end
    
    SNN(expression', 'edge_file.txt', par.k, par.distance);  
    
    system(sprintf('./Cliq.py -i edge_file.txt -o output.txt -r %f -m %f', par.r, par.m));
    SNN_labs = dlmread('output.txt');
    
    unasigned = find(SNN_labs == -1);
    unasigned_labels = max(SNN_labs) + (1:numel(unasigned));
    SNN_labs(SNN_labs == -1) = unasigned_labels;    
end

