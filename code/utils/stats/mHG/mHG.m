function [HGT_vec, pval] = mHG(v, method)
    v = v(:);
	[HGT_vec] = compute_HGT(v);
    min_p = min(HGT_vec);

    N = length(v);
    K = sum(v~=0);
    switch(method)
        case 'exact'
            mat = zeros(K,N-K);	
            pval = mHG_pvalue(N,K,min_p,mat);
        case 'bonferroni'
            pval = min_p*N;
    end
    pval = min(pval, 1);
end

function [pval] = mHG_pvalue(N,K,mHGT,mat)

	if mHGT >= 1
		pval = 1;
		return
	elseif K == 0 || K >= N
		pval = 0;
		return
	end
	
	W = N-K;
	R_ZONE = 0;
	baseHG = 1;
	mat(1,1) = 1;
	
	for n=1:(N-1)
		
		if K >= n
			min_nK = n;
			baseHG = baseHG * (K-n+1) / (N-n+1);
		else
			min_nK = K;
			baseHG = baseHG * n / (n-K);
		end
		
		tailHG = baseHG;
		currHG = baseHG;
		k = min_nK;
		w = n - min_nK;
		
		while tailHG <= mHGT+1e-16 && k > 0 && w < W
			
% 			[n,tailHG,mHGT]
			
			mat(k+1,w+1) = R_ZONE;
			currHG = currHG * (k*(N-K-n+k)) / ((n-k+1)*(K-k+1));
			tailHG = tailHG + currHG;
			
			w = w+1;
			k = k-1;
			
		end
			
		while k >= 0 && w <= W
			
			if w > 0 && k > 0
				mat(k+1,w+1) = mat(k+1,w)*(W-w+1)/(N-n+1) + mat(k,w+1)*(K-k+1)/(N-n+1);
			elseif w > 0
				mat(k+1,w+1) = mat(k+1,w)*(W-w+1)/(N-n+1);
			elseif k > 0
				mat(k+1,w+1) = mat(k,w+1)*(K-k+1)/(N-n+1);
			end
			
			w = w+1;
			k = k-1;
			
		end
		
	end
	pval = 1 - (mat(K+1,W) + mat(K,W+1));
end

function [tail] = HGT(currHG,n,N,K,k)

	min_nK = min(n,K);
	tail = currHG;
	for i=k:(min_nK-1)
		currHG = currHG*((n-1)*(K-i))/((i+1)*(N-n-K+i+1));
		tail = tail+currHG;
	end

end

function [HGT_vec] = compute_HGT(v)

	N = length(v);
	K = sum(v~=0);
    HGT_vec = zeros(N, 1);
    
	currHG = 1;
	k = 0;
		
    for n = 0:N-1
        if v(n+1) == 0
			currHG = currHG*((n+1)*(N-K-n+k))/((N-n)*(n-k+1));
			currHGT = HGT(currHG,n+1,N,K,k);
		else
			currHG = currHG*((n+1)*(K-k))/((N-n)*(k+1));
			k = k+1;
			currHGT = HGT(currHG,n+1,N,K,k);			
        end
        HGT_vec(n+1) = currHGT;
    end
end
