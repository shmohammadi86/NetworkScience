function [hFig] = plotKernel_HeatMap(K, Labels)
    K(K < 0) = 0;
%     K = (K - min(K(:))) ./ (max(K(:)) - min(K(:)));
%     K = 2*K - 1;
    
    %Labels = cellstr(num2str(cell2mat(Labels)));
    n = size(K,1);
    vec = randperm(n);
    K = K(vec,vec);
    Labels = Labels(vec);
    [uniq_Labels] = unique(Labels);
    k = numel(uniq_Labels);
    [~, labels] = ismember(Labels, uniq_Labels);

    cluster_memberships = arrayfun(@(i) sum(labels == i), 1:k);
    [~, cluster_index] = sort(cluster_memberships);
    K_index = arrayfun(@(i) find(labels == cluster_index(i))  ,1:k, 'UniformOutput', false);      

    for i = 1:k
        curr_samples = K_index{i};
        subK = K(curr_samples, curr_samples);
        
        labels = LPA_disjoint(subK);
        [~, order] = sort(labels);
        
        K_index{i} = curr_samples(order);
    end
    K_sort_index = cell2mat(K_index');

    hFig = figure;
    set(hFig, 'Position', [0 0 1000 1000])
    heatmap(K(K_sort_index, K_sort_index));        
    hold all
    
    XLim = get(gca, 'XLim');
    YLim = get(gca, 'YLim');
    num_covered_tissues = cumsum(cluster_memberships(cluster_index))';
    SepLines = num_covered_tissues(1:end-1)';
    alpha = 0.75;
    x1 = repmat([0 XLim(2)]',[1 numel(SepLines)]);
    y1 = repmat(  SepLines(1:end)+YLim(1), [2 1]);
    x2 = repmat(  SepLines(1:end)+XLim(1), [2 1]);
    y2 = repmat([0 YLim(2)]',[1 numel(SepLines)]);
    plot( x1, y1, x2, y2, '-', 'color', [0,0,0]+alpha);
    
    % Setting the Labels for each cluster
    TickLabels = uniq_Labels(cluster_index);
    Ticks = num_covered_tissues(1:k) - (cluster_memberships(cluster_index)' ./ 2);
    set(gca,'YTick', Ticks+YLim(1))
    set(gca, 'YTickLabel', TickLabels, 'Fontsize',22, 'FontWeight','bold');
    set(gca,'fontsize',10);
    set(gca, 'Xtick', Ticks+XLim(1));
    set(gca, 'Xticklabel', TickLabels);
    xticklabel_rotate(Ticks+XLim(1), 60, TickLabels,  'Fontsize',14, 'FontWeight','bold');
%     title(plot_title, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Interpreter', 'None');    
    set(gca,'FontSize', 14, 'FontWeight','bold', 'FontWeight','bold');
    colormap('hot');   
    hold off;
end
