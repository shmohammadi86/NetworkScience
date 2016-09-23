function [ Precision, Recall, FPR, AUC, Max_F1, Max_F1_prec, Max_F1_recall] = validate( Scores, labels, varargin )
%validate: Computes statistics for ranked list of genes coming out of each method.
% Input) Scores: matrix of size number_of_genes*number_of_methods, where number_of_methods is the number of genes and number_of_genes
% is the number of methods, labels: a vector of size number_of_methods*1 where labels(i) =
% 1, if gene i is marked as positive, and 0, otherwise.
% Output) HGT: Hypergeometric tail (Scores) for different cut-sizes., Precision:
% TP/PP = TP/TP+FP. Recall: Also known as TPR/sensitivity/hit rate, is
% defined as TP/P = TP / TP + FN.  FPR: False-Positive Rate, or 1-specificity, defined as
% FP/N = FP/FP+TN. mHG_pval: p-value of the mHG score. Max_F1: max F1 score of the PR curve, where F1 is
% defined as  2*recall*precision/(recall+precision) = 2TP / (2TP + FP +
% FN). AUC: Area-Under-the-Curve of ROC curve.


    params = inputParser;    
    params.addParamValue('correction', 'exact'      ,@(x) ischar(x)); % correction method to use for multiple hypothesis testing
    params.addParamValue('name', 'eval'      ,@(x) ischar(x)); % name of the experiment
    params.addParamValue('method_names', cellstr(arrayfun(@(x) sprintf('Method %d', x), 1:size(Scores, 2), 'UniformOutput', false))      ,@(x) iscell(x) & numel(x) == size(Scores, 2)); % name of each method being evaluated
    params.addParamValue('output_path', './output'      ,@(x) ischar(x)); % output folder to store figures    
    params.parse(varargin{:});
    par = params.Results;
    if(~isvector(par.method_names))
        error('Method names should be a 1D cell array.');
    else
        if(size(par.method_names, 2) == 1)
            par.method_names = par.method_names';
        end
    end
    
    [number_of_genes, number_of_methods] = size(Scores);
     labels = labels(:);
     if(numel(labels) ~= number_of_genes)
        error('number of rows should match in Scores and labels');         
     end

    num_pos = sum(labels);
    if num_pos==0
        error('no positive labels entered');
    end
    if num_pos==number_of_genes
        error('no negative labels entered');
    end
    
    if (~exist(par.output_path, 'dir'))
        mkdir(par.output_path);
    end
    
    perm = randperm(size(Scores, 1), size(Scores, 1));
    Scores = cell2mat(arrayfun(@(c) Scores(perm, c), 1:size(Scores, 2), 'UniformOutput', false));
    labels = labels(perm);

    [~, scores_si] = sort(Scores,'descend');
    %clear Scores
    scores_si_reindex = scores_si+ones(number_of_genes,1)*(0:number_of_methods-1)*number_of_genes;
    Labels = repmat(labels, 1, number_of_methods);
    l = Labels(scores_si_reindex);
    clear scores_si Labels 

    tp = cumsum(l==1,1);
    fp = repmat((1:number_of_genes)',[1 number_of_methods])-tp;

    num_neg = number_of_genes-num_pos;
    Precision = bsxfun(@rdivide,tp, (1:number_of_genes)'); 
    Recall = bsxfun(@rdivide,tp,num_pos); 
    FPR = bsxfun(@rdivide,fp,num_neg); 

    F1 = 2*Recall .* Precision ./ (Recall + Precision);
    [Max_F1, Max_F1_idx] = max(F1);
    Max_F1_prec = Precision(Max_F1_idx);
    Max_F1_recall = Recall(Max_F1_idx);
    
    AUC = sum(Recall.*[(diff(fp)==1); zeros(1, number_of_methods)])./num_neg;

%     % Plot PR curve
%     fig = figure;
%     set(fig, 'PaperOrientation', 'portrait');
%     set(fig,'units','normalized','outerposition',[0 0 1 1]);
%     
%     plot([zeros(1, number_of_methods); Recall], [ones(1, number_of_methods); Precision], 'LineWidth', 2);        
%     set(gca,'FontSize', 16, 'FontWeight','bold'); 
%     xlabel('Recall', 'FontSize',  17, 'FontWeight', 'bold'); 
%     ylabel('Precision', 'FontSize', 17, 'FontWeight', 'bold');
% 
%     legend(cellfun(@(x, y) sprintf('%s (F1 = %.3f)', x, y), par.method_names, num2cell(full(Max_F1)), 'UniformOutput', false), 'FontSize', 18, 'FontWeight', 'bold', 'interpreter', 'none');
%     xlim([0 1]); ylim([0 1]);    
%     
%     full_path = fullfile(par.output_path, sprintf('%s_PR.pdf',  par.name));    
%     export_fig(full_path, '-pdf', '-transparent');
%     close;    

    % Plot ROC curve
    fig = figure;
    set(fig, 'PaperOrientation', 'portrait');
    set(fig,'units','normalized','outerposition',[0 0 1 1]);
    
    plot([zeros(1, number_of_methods); FPR; ones(1, number_of_methods)], [zeros(1, number_of_methods); Recall; ones(1, number_of_methods)], 'LineWidth', 2);        
    hold on;
    plot(0:0.01:1, 0:0.01:1, '--k')    
    hold off;
    set(gca,'FontSize', 16, 'FontWeight','bold'); 
    xlabel('1 - Specificity', 'FontSize', 17, 'FontWeight', 'bold'); 
    ylabel('Sensitivity', 'FontSize', 17, 'FontWeight', 'bold');
    legend(cellfun(@(x, y) sprintf('%s (AUC = %.3f)', x, y), par.method_names, num2cell(full(AUC)), 'UniformOutput', false), 'FontSize', 18, 'FontWeight', 'bold', 'interpreter', 'none', 'Location', 'SouthEast');
    
    xlim([0 1]); ylim([0 1]);    
    
    title(par.name);
    full_path = fullfile(par.output_path, sprintf('%s_ROC.pdf',  par.name));    
    export_fig(full_path, '-pdf', '-transparent');
    close;    
    
    
end
