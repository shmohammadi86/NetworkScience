function [ raw_transcriptional_signatures, raw_transcriptional_similarity, normalized_common_signature, adjusted_transcriptional_signatures,  adjusted_transcriptional_similarity ] = peelMat( T, varargin )
%peelMat 
    params = inputParser;    
    params.addParamValue('type', 'SVD'      ,@(x) ischar(x)); 
%     params.addParamValue('transform', true      ,@(x) islogical(x)); 
    params.parse(varargin{:});
    par = params.Results;
    
    probeset_count = size(T, 1);
    raw_transcriptional_signatures = zscore(T); 
    raw_transcriptional_similarity =  raw_transcriptional_signatures' * raw_transcriptional_signatures ./ probeset_count; %Eq: corr(T); i.e., raw transcriptional similarity
    
    switch(par.type)
        case 'mean'
            common_signature = mean(raw_transcriptional_signatures, 2) ; 
            

        case 'max'
            common_signature = max(raw_transcriptional_signatures, [], 2);             
        case 'SVD'
            [U,S,V] = svds(raw_transcriptional_signatures, 1); 
            if sum(V(:, 1)) < 0
                U = -U; 
                V = -V;
            end
            sigma = diag(S); 
            common_signature = sigma(1) * U(:, 1);
    end                
    normalized_common_signature = zscore(common_signature); % First-order approximation of common (i.e. housekeeping gene) signatures
    
    
    common_correlation = normalized_common_signature'*raw_transcriptional_signatures ./ probeset_count;  %Eq: corr(T, normalized_common_signature)'; i.e., correlation between each celltype and the common signature
    adjusted_T  = raw_transcriptional_signatures - normalized_common_signature*common_correlation;    
    adjusted_transcriptional_signatures = zscore(adjusted_T); %adjusted_transcriptional_signatures 

    adjusted_transcriptional_similarity = adjusted_transcriptional_signatures'*adjusted_transcriptional_signatures ./ probeset_count; % Eq:  partialcorr(T, common_signature); i.e., adjusted  transcriptional similarity 

%     if(par.transform)
%         raw_transcriptional_signatures = Sigmoid(raw_transcriptional_signatures);
%         adjusted_transcriptional_signatures = Sigmoid(adjusted_transcriptional_signatures);
%         normalized_common_signature = Sigmoid(normalized_common_signature);
%     end
%     switch(type)
%         case 'mean'
%             common_signature = mean(T, 2); 
%         case 'SVD'
%             % Check compatibility
%             residual_T = T - (U(:, 1)*S(1, 1)*V(:, 1)');            
%             corr(mean(T, 2), common_signature, 'type', 'spearman')
%             residual_corr = corr(residual_T);
%             corr(residual_corr(:), adjusted_transcriptional_similarity(:))
%     end                    
end

