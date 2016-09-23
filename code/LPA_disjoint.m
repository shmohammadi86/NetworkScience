function [ labels ] = LPA_disjoint( A, varargin )
%% A) Similarity matrix (adjacency matrix of a graph) of data points
% Returns cell array of clusters

    [~, n] = size(A); % gene * sample: n = #samples    
    params = inputParser;    
    params.addParamValue('max_it',  30, @(x) isscalar(x) & x > 0);
    params.addParamValue('rounds', 10, @(x) isscalar(x) & x > 0);
    params.addParamValue('ratio', 1, @(x) isscalar(x) & (0 <= x & x <= 1) | x == -1 ); % non-overlapping    
    params.addParamValue('labels', {} , @(x) iscell(x) & (numel(x) == n | isempty(x)) );
    params.addParamValue('verbose', false, @(x) islogical(x));
    params.parse(varargin{:});
    par = params.Results;
    
    A(A < 0) = 0;
    A = A - diag(diag(A));
    if(nnz(A) == 0)
        labels = 1:size(A, 1);
        return;
    end    
    A = A ./ max(nonzeros(A));
    
    deg = sum(A);
    
    standard_labels = cell(par.rounds, 1);    
    tail = n*ones(5, 1); tail_idx = 0;
    for i = 1:par.rounds
        if(par.verbose)
            fprintf('Round %d\n', i);
        end
        
        % Initialize node clusters
        k = randsample(round(n/3) ,1); % do we need this?
        l = randsample(k, n, true);
        L = sparse(1:n, l, 1, n, k); % samples * clusters
        deg_rep = repmat(deg, k, 1)';


        % Main label propagation loop
        for it = 1:par.max_it
            if(par.verbose)
                fprintf('\titeration %d\n', it);
            end

            % Adjust for global frequency of labels to allow more refined clusters (instead of one label over-reaching all over the network)
            global_freq = mean(L);
            expected_local_freq = sparse(bsxfun(@times, deg_rep, global_freq)); %expected_local_freq = deg'*global_freq;
            
            L_prime = A*L - expected_local_freq;
            
            % Update L (labels matrix) for the next round
            [~, new_l] = max(L_prime, [], 2);            
            L = sparse(1:numel(l), new_l, 1, n, k);                             
            
            % For faster computation
            changes = nnz(l ~= new_l);
            tail(tail_idx+1) = changes;
            tail_idx = mod(tail_idx+1, 5);
            changes_mean = mean(tail);
            if(changes_mean < round(0.01*n))
                break;
            end            
            l = new_l;
        end
        
        % Convert cluster ids to standard form 1 ... k
        [~, standard_labels{i}] = ismember(l, unique(l));

    end

    overlap = zeros(par.rounds);
    for i = 1:par.rounds
        for j = i+1:par.rounds
            overlap(i, j) = adjustedrand(standard_labels{i}, standard_labels{j});
        end
    end
    overlap = max(overlap, overlap');
    [~, core_cluster_id] = max(sum(overlap));
    
    labels = standard_labels{core_cluster_id};    
end

