function TF = idx2logical(idxs,n)
    if all(ismember(idxs,[0,1]))
        mustBeVector(idxs)
        mustBeLessThanOrEqual(numel(idxs),n)
        TF = logical(idxs(:));
    else
        mustBeLessThanOrEqual(idxs,n)
        TF = false(n,1);
        TF(idxs) = true;
    end
end