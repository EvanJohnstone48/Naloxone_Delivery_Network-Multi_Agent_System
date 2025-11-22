function centroids = compute_coverage_centroids(P, X, Y, D)

    if isempty(P)
        centroids = [];
        return;
    end

    [nRows, nCols] = size(D);
    nCov  = size(P, 1);

    xVec = X(:);
    yVec = Y(:);
    wVec = D(:);

    nCells = numel(wVec);
    owner  = zeros(nCells, 1);

    % Assign each cell to nearest coverage drone
    for k = 1:nCells
        dx = P(:,1) - xVec(k);
        dy = P(:,2) - yVec(k);
        distSq = dx.^2 + dy.^2;
        [~, idxMin] = min(distSq);
        owner(k) = idxMin;
    end

    % Compute centroids
    centroids = zeros(nCov, 2);
    for j = 1:nCov
        mask = (owner == j);
        wj   = wVec(mask);

        if isempty(wj) || sum(wj) == 0
            centroids(j,:) = P(j,:);
        else
            xj = xVec(mask);
            yj = yVec(mask);
            wSum = sum(wj);
            cx = sum(xj .* wj) / wSum;
            cy = sum(yj .* wj) / wSum;
            centroids(j,:) = [cx, cy];
        end
    end
end
