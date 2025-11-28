function inspect_toronto_geojson()
    % Simple inspector for your toronto.geojson -> shapes struct

    clc;
    fprintf('--- Inspecting toronto.geojson ---\n');

    shapes = load_toronto_geojson();

    nShapes = numel(shapes);
    fprintf('Number of shapes (entries): %d\n', nShapes);

    if nShapes == 0
        fprintf('shapes is EMPTY. That is the problem.\n');
        return;
    end

    % We built load_toronto_geojson as a single point cloud:
    allLon = shapes.lon(:);
    allLat = shapes.lat(:);

    fprintf('Raw lon range: [%.6f, %.6f]\n', min(allLon), max(allLon));
    fprintf('Raw lat range: [%.6f, %.6f]\n', min(allLat), max(allLat));

    % Keep only values that look like Toronto:
    %   lon ~ [-80, -78], lat ~ [43, 45].
    validMask = (allLon > -80 & allLon < -78) & ...
                (allLat >  43 & allLat <  45);

    lonT = allLon(validMask);
    latT = allLat(validMask);

    fprintf('Filtered to %d plausible Toronto points.\n', numel(lonT));
    if isempty(lonT)
        warning('No points passed the Toronto filter. Check coordinate ranges.');
    else
        fprintf('Filtered lon range: [%.6f, %.6f]\n', min(lonT), max(lonT));
        fprintf('Filtered lat range: [%.6f, %.6f]\n', min(latT), max(latT));
    end

    figure(123); clf;
    if isempty(lonT)
        text(0.5, 0.5, 'No filtered points to show', ...
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle', ...
             'FontSize', 14);
        axis off;
        return;
    end

    plot(lonT, latT, 'r.', 'MarkerSize', 6);
    xlabel('Longitude'); ylabel('Latitude');
    title('Toronto neighbourhood point cloud (filtered)');
    axis equal;
    xlim([-80 -78]);
    ylim([43 45]);
    grid on;
end

