function [ Graph ] = my_readGraph( fname )
    fd = fopen(fname, 'r');
    stats = textscan(fd, '%f %f', 1, 'Whitespace', '', 'Delimiter', '\t');
    n = stats{1};
%     edges = textscan(fd, '%f %f %.16f', 'Whitespace', '', 'Delimiter', '\t');
    edges = fscanf(fd, '%f\t%f\t%f', [3, Inf])';
    fclose(fd);
       
    Graph = sparse(edges(:, 1), edges(:, 2), edges(:, 3), n, n);
    Graph = max(Graph, Graph') - diag(diag(Graph)); % Remove self-loops, if exists, and convert to undirected graph (symmetrize)
end

