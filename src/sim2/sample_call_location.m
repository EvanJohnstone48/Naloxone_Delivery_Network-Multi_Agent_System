function [cx, cy, iRow, iCol] = sample_call_location(X, Y, D)
    w = D(:);
    w = w - min(w); 
    if sum(w) <= 0
        error('sample_call_location:ZeroDensity', ...
              'Density map has no positive entries.');
    end
    w = w / sum(w);          % normalize to probabilities

    cdf = cumsum(w);
    r = rand;
    idx = find(cdf >= r, 1, 'first');

    cx = X(idx);
    cy = Y(idx);
    [iRow, iCol] = ind2sub(size(D), idx);
end
