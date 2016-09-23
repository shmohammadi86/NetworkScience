function [label, energy] = k2means(K, init)
n = size(K,2);
if numel(init)==1
    k = init;
    label = ceil(k*rand(1,n));
elseif numel(init)==n
    label = init;
end


last = 0;
max_iter = 100;
it = 0;
while any(label ~= last) && it < max_iter
    it = it + 1;
    [u,~,label(:)] = unique(label);   % remove empty clusters
    k = numel(u);
    E = sparse(label,1:n,1,k,n,n);
    E = spdiags(1./sum(E,2),0,k,k)*E;
    T = E*K;
    last = label;
    [val, label] = max(bsxfun(@minus,T,diag(T*E')/2),[],1);
end
energy = trace(K)-2*sum(val); 

