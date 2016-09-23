function spyc(sA,cmap)

%SPYC Visualize sparsity pattern with color-coded scale.
%   SPYC(S) plots the color-coded sparsity pattern of the matrix S.
%
%   SPYC(S,CMAP) plots the sparsity pattern of the matrix S USING 
%                    COLORMAP CMAP.
%
%   SPYC(S,CMAP,PB) allows turning off the display of a colorbar by passing
%                   flag PB=0
%
%   written by Try Hard
%   $Revision: 0.0.0.2 $  $Date: 2013/08/24 11:11:11 $

if nargin<1 | nargin>3 && ~isempty(cmap)
    error( 'spyc:InvalidNumArg', 'spyc takes one to three inputs')
    return
end

if isempty(sA)
    error( 'spyc:InvalidArg', 'sparse matrix is empty')
    return
end

if nargin>1 && ~isempty(cmap)
    % colorspy does not check whether your colormap is valid!
    if ~isnumeric(cmap)
        cmap=flipud(colormap(cmap));
    end
else
    cmap=flipud(colormap('gray'));
end


sA = sA ./ sum(sA(~isnan(sA) & ~isinf(sA)));
[Nx Ny]=size(sA);

NNZ = round(sum(sA(~isnan(sA) & ~isinf(sA)))^2/sum(sA(~isnan(sA) & ~isinf(sA)).^2));
%NNZ = nnz(sA > mean(sA(~isnan(sA) & ~isinf(sA))));
% sA(sA < mean(sA(:))) = 0;
% if(nnz(sA) > 1e4)
    [~, sA_idx] = sort(sA(:), 'descend');
    sA(sA_idx(NNZ:end)) = 0;
% end
indx=find(sA);
sA=full(sA(indx));
ns = length(indx);
[ix iy]=ind2sub([Nx Ny],indx);

if(max(sA)-min(nonzeros(sA)) == 0)
    imap = 63*ones(size(sA));
else
    imap = round((sA-min(nonzeros(sA)))*63/(max(sA)-min(nonzeros(sA))))+1;
end

% figure, hold on
colormap(cmap)
scatter(iy,ix,[],imap,'Marker','.','SizeData',200)
set(gca,'ydir','reverse')
axis equal;
axis([0 Ny 0 Nx])
box on

set(gcf,'units','normalized','outerposition',[0 0 1 1]); 



