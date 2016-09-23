function [mi,mj, val] = matching_sparse(L)
% MATCHING_SPARSE
%
% [mi,mj,val] = matching_sparse(X)
%
% Solve the bipartite, maximum weight matching problem on the non-zero
% entries of L and return the matching subset.

	[n m] = size(L);
	[i j v] = find(L);
	nedges = nnz(L);

	[mi mj matchsize val] = matching_sparse_mex(n,m,nedges,i-1,j-1,v);

	mi = mi(1:matchsize) + 1;
	mj = mj(1:matchsize) + 1;
end

