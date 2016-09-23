clear 

addpath(genpath('./code/'));

paths.output = 'output/part1/';
if(~exist(paths.output, 'dir'))
    mkdir(paths.output);
end

ds_path = 'input/transcriptome/Pollen';


%% Read expression dataset
    tic; [expression, sample_names, gene_names] = my_tblread(fullfile(ds_path, 'expression.txt')); toc
    unexpressed_mask = sum(expression, 2) == 0;
    expression(unexpressed_mask, :) = [];
    gene_names(unexpressed_mask) = [];

%% Import cell labels
    sample_annotations = my_dlmread(fullfile(ds_path, 'sample_annotations.txt'));

    Labels = sample_annotations(:, end);
    k = numel(unique(Labels));
    UL = unique(Labels);
    [~, true_labels] = ismember(Labels, UL);
    
  
%% Compute cell similarities
    [ raw_transcriptional_signatures, raw_transcriptional_similarity, normalized_common_signature, adjusted_transcriptional_signatures,  adjusted_transcriptional_similarity ] = peelMat( expression );       
    
    
%% Plot heat map of similarity matrix before/after adjustment
    
    plotKernel_HeatMap(raw_transcriptional_similarity, Labels);

    plotKernel_HeatMap(adjusted_transcriptional_similarity, Labels);

    
%% Plot 2D projection of points
    rng('default')
    
    ydata = tsne_d(1-corr(adjusted_transcriptional_signatures), Labels, 2);   
    

    fig = figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1], 'PaperPositionMode', 'auto');  
    CC = cbrewer('qual', 'Set1', numel(UL));
    hold all
    for i = 1:numel(UL)
        idx = ismember(Labels, UL{i});
        scatter(ydata(idx,1), ydata(idx,2), 25, CC(i, :), 'filled');
    end
    hold off
    legend(UL);   
    set(gca,'FontSize', 16, 'FontWeight','bold'); 
    xlabel('x-tsne');
    ylabel('y-tsne');
    

%% Distribution of similrity scores
    figure
    ksdensity(adjusted_transcriptional_similarity(:))
    xlabel('Similarity');
    ylabel('pdf');
    
%% Construct tissue-tissue similarity network (TTSN)
    [ii, jj] = find(tril(adjusted_transcriptional_similarity, -1) > 0.4);
    
    vv = adjusted_transcriptional_similarity(sub2ind(size(adjusted_transcriptional_similarity), ii, jj));
    
    Edgelist = [sample_annotations(ii, 1), sample_annotations(jj, 1), num2cell(vv)];
    dlmcell(fullfile(paths.output, 'TTSN_edgeList.txt'), Edgelist);
    
    node_annotations = [sample_annotations(:, [1, 3, 5])];
    NodeAnnotations = [{'ID', 'Cell', 'Tissue'}; node_annotations];
    dlmcell(fullfile(paths.output, 'TTSN_nodeAnnotations.txt'), NodeAnnotations);
    
%% Go to Cytoscape and visualize the TTSN

    
%% Validate cell type identification

    % kmeans
    sample_no = 100; 
    NMI_samples = zeros(sample_no, 1);
    ARI_samples = zeros(sample_no, 1);
    rng('default');            
    for j = 1:sample_no
        fprintf('\tsample %d/%d\n', j, sample_no);
        labels = k2means(adjusted_transcriptional_similarity, k);
        NMI_samples(j) = nmi(labels, true_labels);
        ARI_samples(j) = adjustedrand(labels, true_labels);
    end

    mean_NMI = mean(NMI_samples);
    mean_ARI = mean(ARI_samples);            
    
    % LPA
    LPA_labels = LPA_disjoint(adjusted_transcriptional_similarity);
    LPA_NMI = nmi(LPA_labels, true_labels);
    LPA_ARI = adjustedrand(LPA_labels, true_labels); 
    
    
    % SNN_Cliq
    cd('code/SNNCliq/');
    SNN_labels = run_SNNCliq( expression );
    cd('../..');
    
    SNN_NMI = nmi(SNN_labels, true_labels);
    SNN_ARI = adjustedrand(SNN_labels, true_labels);     
    

%%


    
%% Identify cell type-specific markers
    
    curr_cluster = 4;
    neural_mask = (LPA_labels == curr_cluster);

    pvals = mattest(expression(:, neural_mask), expression(:, ~neural_mask));    
    [~,~,padj] = fdr(pvals);
    
    fold = log2(mean(expression(:, neural_mask), 2) ./ mean(expression(:, ~neural_mask), 2));
    neural_genes = gene_names(padj < 1e-3 & fold > 2);
    

    % compute marker genes
    [~, perm] = sort(LPA_labels);    
    neural_mask_perm = neural_mask(perm);
    
    unique_clusters = unique(LPA_labels);
    NRREPLICS = arrayfun(@(idx) nnz(LPA_labels == unique_clusters(idx)), 1:numel(unique_clusters));
    
   
    [BERGER, IUT_PVAL]=IUT(expression(:, perm),NRREPLICS,curr_cluster,1,0.05,1);
    [~,~,iut_padj] = fdr(IUT_PVAL);    
    neural_genes_IUT = gene_names(iut_padj < 1e-3 & fold > 2);

    
    % Bayes factor
    r=Constr_Bayes(expression(:, perm),NRREPLICS,curr_cluster,1, 1000);
    neural_genes_Bayes = gene_names(r > 100);
      
%% Do Enrichment analysis on markers

    