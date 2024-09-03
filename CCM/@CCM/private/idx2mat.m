function TF = idx2mat(idxs,n)
    if all(ismember(idxs,[0,1]))
        mustBeVector(idxs)
        mustBeLessThanOrEqual(numel(idxs),n)
        idxs = find(idxs == 1);
    else
        mustBeLessThanOrEqual(idxs,n)
    end
    m = numel(idxs);
    TF = false(m,n);
    idxs = sub2ind(size(TF), 1:m, reshape(idxs(:),1,m));
    TF(idxs) = true;
end