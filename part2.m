clear 

%%
addpath(genpath('./code/'));


paths.output = 'output/part2/';
if(~exist(paths.output, 'dir'))
    mkdir(paths.output);
end

paths.transcriptome = 'input/transcriptome/Pollen';
paths.interactome = 'input/interactome';


%% Read transcriptome
    tic; [expression, sample_names, gene_names] = my_tblread(fullfile(paths.transcriptome, 'expression.txt')); toc
    unexpressed_mask = sum(expression, 2) == 0;
    expression(unexpressed_mask, :) = [];
    gene_names(unexpressed_mask) = [];

    sample_annotations = my_dlmread(fullfile(paths.transcriptome, 'sample_annotations.txt'));

    
%% Load interactome and filter to transcriptome genes
    
    node_annotations = readtable(fullfile(paths.interactome, 'human_interactome_node_annotations.txt'), 'Delimiter', '\t');

    
    common_genes= intersect(unique(node_annotations.hgnc_symbol), gene_names);    
    
    expression(~ismember(gene_names, common_genes), :) = [];
    gene_names(~ismember(gene_names, common_genes)) = [];
    [~, perm] = sort(gene_names);
    gene_names = gene_names(perm);
    expression = expression(perm, :);
    

    [~, idx] = ismember(common_genes, node_annotations.hgnc_symbol);
    node_annotations = node_annotations(idx, :);
    
    

    interactome_edge_list = readtable(fullfile(paths.interactome, 'human_interactome_edge_list_simplified.txt'), 'Delimiter', '\t');        
    [ii_mask, ii] = ismember(interactome_edge_list.icrogidA, node_annotations.icrogid);
    [jj_mask, jj] = ismember(interactome_edge_list.icrogidB, node_annotations.icrogid);    
    
    node_mask = ii_mask & jj_mask;
    
    interactome = sparse(ii(node_mask), jj(node_mask), interactome_edge_list.w(node_mask), size(node_annotations, 1), size(node_annotations, 1), nnz(node_mask));
    interactome = max(interactome, interactome'); % symmetrize       


%% Indetify blood-specific pathways/complexes
%     [~, perm] = sort(sum(interactome), 'descend');
%     party_genes = gene_names(perm(1:30))
    unspecific_hubs = ismember(gene_names, {'UBC', 'SUMO1', 'SUMO2', 'TP53'});
    pruned_interactome = interactome;
    pruned_interactome(unspecific_hubs, :) = 0;
    pruned_interactome(:, unspecific_hubs) = 0;
    
    mc_samples = 1000;
    r_threshold = 10;
    
    fprintf('Compute Bayes Factors with %d samples\n', mc_samples);
    [~, perm] = sort(sample_annotations(:, 5));    
    
    ut = unique(sample_annotations(:, 5));
    NRREPLICS = arrayfun(@(idx) nnz(strcmp(ut{idx}, sample_annotations(:, 5))), 1:numel(ut));    
    r=Constr_Bayes(expression(:, perm),NRREPLICS,find(strcmp(ut, 'BLOOD')),1, mc_samples);    

    blood_markers = gene_names(r > r_threshold);

    fprintf('Construct blood-specific interactome\n');
    blood_adj_net = construct_TS_interactome(pruned_interactome, r, 'method', 'PenPro_L1', 'expression_threshold', r_threshold, 'lambda', 0.1);
    blood_adj_net(blood_adj_net < 1) = 0;
    
    
    [ii, jj] = find(tril(blood_adj_net, -1));
    vv = blood_adj_net(sub2ind(size(blood_adj_net), ii, jj));
    vv = normalize(log(vv), 'pnorm', 1)*numel(vv);
    
    
    blood_edge_list = [gene_names(ii), gene_names(jj), num2cell(vv)];
    
    dlmcell('output/part2/blood_diffNet_edgeList.txt', blood_edge_list);

%% Load in Cytoscape and do clustering using MCL, then do enrichment of clusters


%% Loading markers    
    fprintf('Loading leukemia marker genes \n');
    % from https://www.intogen.org/search?cancer=AML
    Leukemia_genes = {'DNMT3A', 'FLT3', 'RUNX1', 'TP53', 'IDH2', 'IDH1', 'NRAS', 'KRAS', 'U2AF1', 'NPM1', 'TET2', 'CEBPA', 'WT1', 'KIT', 'PTPN11', 'STAG2', 'PHF6', 'RAD21', 'ASXL1', 'KDM6A', 'CHD4', 'EGFR', 'NF1', 'EZH2', 'DIS3', 'SUZ12', 'THRAP3', 'CUL1', 'CBFB', 'PRPF8', 'BCOR', 'MED12'};
    
    [mask, gene_idx] = ismember(Leukemia_genes, common_genes);
    Leukemia_genes(~mask) = [];
    gene_idx(~mask) = [];
    
    gold = sparse(gene_idx, 1, 1, numel(common_genes), 1);


%% Preprocess
    fprintf('Normalize scores between [0, 1]\n');
    sigmoid_raw = 1 ./ (1 + exp(-zscore(expression)));

    leukemia_mask = strcmp(sample_annotations(:, 3), 'K562');
    leukemia_submat = sigmoid_raw(:, leukemia_mask);
    leukemia_signature = max(leukemia_submat, [], 2);

    
    leukemia_w = normalize(leukemia_signature, 'pnorm', 1)*numel(leukemia_signature);
    ksdensity(leukemia_w)
    
%%  Run disease gene prioritization
    pvals = ones(numel(gene_idx), 2);    

    fprintf('Running disease gene prioritization\n');
    
    fprintf('\Running on global network\n');
    for i = 1:numel(gene_idx)
        fprintf('\t\tGene %d/%d\n', i, numel(gene_idx));
        src_gene_idx = gene_idx(i);
        v = sparse(src_gene_idx, 1, 1, numel(common_genes), 1);
        
        pr = pagerank(interactome, struct('v', v, ...
                                                  'alg','linsys',... 
                                                  'linsys_solver',@(f,v,tol,its) gmres(f,v,[],tol, its)));                       
        [~, perm] = sort(pr, 'descend');
        [~, pvals(i, 1)] = mHG(gold(perm), 'exact');
    end

           

    fprintf('\tConstructing leukemia-specific network\n');
    leukemia_net = construct_TS_interactome(interactome, leukemia_w, 'method', 'PenPro_L1', 'expression_threshold', 0.8, 'alpha', 0.15);
    fprintf('\tRunning on tissue-specific network\n');
    for i = 1:numel(gene_idx)
        fprintf('\t\tGene %d/%d\n', i, numel(gene_idx));
        src_gene_idx = gene_idx(i);
        v = sparse(src_gene_idx, 1, 1, numel(common_genes), 1);
        
        pr = pagerank(leukemia_net, struct('v', v, ...
                                                  'alg','linsys',... 
                                                  'linsys_solver',@(f,v,tol,its) gmres(f,v,[],tol, its)));                       
        [~, perm] = sort(pr, 'descend');
        [~, pvals(i, 2)] = mHG(gold(perm), 'exact');
    end
    
    combined_pvals = arrayfun(@(col) Edgington(pvals(:, col)), 1:size(pvals, 2))
    

%% Construct Disease-related pathways
    lambda = 1;% Importance of node prizes compared to edge penalties (smaller the value of lambda, smaller the size of constructed pathway)
    depth = 4; % max depth of the rooted PCST 
    max_iter = 10;

    fprintf('Construct Disease-related pathways\n');

    fprintf('Remove party hubs\n');
    party_hubs = find(sum(logical(interactome)) > 450);
    Pruned_net = leukemia_net;
    Pruned_net(party_hubs, :) = 0;
    Pruned_net(:, party_hubs) = 0;   

    fprintf('Construct resistance network \n');
    [ii, jj, vv] = find(tril(Pruned_net));
    transformed_weights = spfun(@(x) 1 ./ x, vv); % Inverse of w_ij / sqrt(deg_i*deg_j)    
    transformed_interactome = sparse(ii, jj, transformed_weights, size(interactome, 1), size(interactome, 1));
    transformed_interactome = max(transformed_interactome, transformed_interactome');    

    
    
    
    % Node prizes!
    fprintf('Compute node prizes \n');
    pruned_leukemia_idx = gene_idx;
    pruned_leukemia_genes = Leukemia_genes;
    
    %remove isolated nodes
    pruned_deg = sum(logical(transformed_interactome(pruned_leukemia_idx, :)), 2);
    pruned_leukemia_idx(pruned_deg == 0) = [];   
    pruned_leukemia_genes(pruned_deg == 0) = []; 

    node_prizes = full( (sum(interactome(pruned_leukemia_idx, :)) ./ sum(interactome)) );
    node_prizes(isnan(node_prizes)) = 0;
    node_prizes(pruned_leukemia_idx) = 1000; % Make sure to Include disease genes


    fprintf('Compute a bunch of trees rooted on different leukemia-related genes\n');
    min_cost = inf;
    for root_id = 1:numel(pruned_leukemia_idx)
        fprintf('Expanding from %s: gene %d/%d\n', pruned_leukemia_genes{root_id}, root_id, numel(pruned_leukemia_genes));
        root = pruned_leukemia_idx(root_id);
        [parents, total_cost] = MsgSteiner(transformed_interactome, lambda*node_prizes, depth, max_iter, root);	
        if(total_cost < min_cost)
            min_cost = total_cost;
            best_tree = parents;
            fprintf('\t\t\tBest: size = %d, cost = %f\n', nnz(parents), total_cost);
        end
    end    

    [~, ii, jj] = find(best_tree);
    idx = sub2ind(size(interactome), max(ii, jj), min(ii, jj));
        
    
    % Visualize in Cytoscape now!
    all_nodes = union(ii, jj);
    node_table = table2cell(node_annotations(all_nodes, 3:end));
    node_table(:, 4) = num2cell(leukemia_w(all_nodes));
    node_table(:, 5) = num2cell(ismember(all_nodes, pruned_leukemia_idx));
    dlmcell('output/part2/Leukemia_Steiner_node_annotations.txt', node_table);

    [local_ii, local_jj, vv] = find(tril(leukemia_net(all_nodes, all_nodes)));
    global_ii = all_nodes(local_ii);
    global_jj = all_nodes(local_jj);
    global_idx = sub2ind(size(interactome), global_ii, global_jj);

    edge_table = [node_annotations.hgnc_symbol(global_ii), node_annotations.hgnc_symbol(global_jj), num2cell(vv), num2cell(zeros(numel(vv), 1))];
    edge_table(ismember(global_idx,idx), 4) = {1};
    dlmcell('output/part2/Leukemia_Steiner_edge_table.txt', edge_table);
        