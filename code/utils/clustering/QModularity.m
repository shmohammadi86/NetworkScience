function [ Q ] = QModularity( A, clusters )
%     A(A < 0) = 0;
    A = A - diag(diag(A));
    
    d_i = sum(A);
    Deg_prod = d_i'*d_i;

    w = sum(d_i);
    scaled_Deg_prod = Deg_prod ./ w;
    
    Delta = A - scaled_Deg_prod;
    Delta = Delta - diag(diag(Delta));
    
    Q = sum(arrayfun(@(id) sum(nonzeros(Delta(clusters==id, clusters==id))), 1:numel(unique(clusters)))) / w;    

end

