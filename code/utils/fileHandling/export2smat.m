function [  ] = export2smat( adj, path)
%% Exports adjacency matrix to smat format
    n = size(adj, 1);
    [ii, jj] = find(tril(adj, -1));
     
     edge_list = [ii-1, jj-1];
    m = numel(ii);
    fd = fopen(path, 'w');
    fprintf(fd, '%d\t%d\t%d\n', n, n, m);
    fprintf(fd, '%d\t%d\t1\n', edge_list');
    fclose(fd);
end

