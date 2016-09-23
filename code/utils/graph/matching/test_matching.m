function test_matching
X = rand(10,3) + 10*eye(10,3);
[mi,mj] = matching_sparse(X);
assert(all(mi == mj));
