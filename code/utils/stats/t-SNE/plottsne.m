function plottsne(S, Labels, plot_title, ds_path)
    unique_labels = unique(Labels);
    [~, labels] = ismember(Labels, unique_labels);        
    ydata = tsne_p(S, labels, 2);
    
    figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1], 'PaperPositionMode', 'auto');  
    set(gcf,'units','normalized','outerposition',[0 0 1 1], 'PaperPositionMode', 'auto');  
    CC = cbrewer('qual', 'Set1', numel(unique_labels));

    hold on;
    for i = 1:numel(unique_labels)
        idx = (labels == i);
        scatter(ydata(idx,1), ydata(idx,2), 25, CC(i, :), 'filled');
    end
    hold off; 

    legend(unique_labels);   
    set(gca,'FontSize', 16, 'FontWeight','bold'); 
    xlabel('x-tsne-pca');
    ylabel('y-tsne-pca');
    title(plot_title, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Interpreter', 'None');
    
    export_fig(fullfile(ds_path, strcat(plot_title, '_tsne.eps')), '-eps', '-transparent', '-painters');
    close;
end