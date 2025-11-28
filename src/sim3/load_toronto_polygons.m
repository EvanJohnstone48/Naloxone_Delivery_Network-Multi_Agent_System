function polys = load_toronto_polygons()
% load_toronto_polygons  Approximate Toronto outline polygon from GeoJSON.
%
% Reads data/toronto_area/toronto.geojson, collects ALL coordinates from
% neighbourhood geometries, and builds a concave boundary polygon.
%
% Output:
%   polys : struct array with fields
%           polys(k).lon, polys(k).lat  (here k=1, whole city outline)

    % ---- Locate toronto.geojson relative to this file ----
    thisFile   = mfilename('fullpath');      % .../src/sim3/load_toronto_polygons.m
    sim3Folder = fileparts(thisFile);        % .../src/sim3
    srcFolder  = fileparts(sim3Folder);      % .../src
    projRoot   = fileparts(srcFolder);       % .../Naloxone...
    dataFolder = fullfile(projRoot, 'data', 'toronto_area');
    geojsonPath = fullfile(dataFolder, 'toronto.geojson');

    if ~isfile(geojsonPath)
        error('load_toronto_polygons:NoFile', ...
              'GeoJSON not found at %s', geojsonPath);
    end

    % ---- Decode JSON ----
    txt = fileread(geojsonPath);
    G   = jsondecode(txt);

    % ---- Collect ALL coordinate pairs into one big cloud ----
    lonAll = [];
    latAll = [];

    if isfield(G, 'features')
        for f = 1:numel(G.features)
            feat = G.features(f);
            if ~isfield(feat, 'geometry') || isempty(feat.geometry)
                continue;
            end
            geom = feat.geometry;
            if ~isfield(geom, 'coordinates')
                continue;
            end
            [lonAll, latAll] = collect_coords(geom.coordinates, lonAll, latAll);
        end
    elseif isfield(G, 'coordinates')
        [lonAll, latAll] = collect_coords(G.coordinates, lonAll, latAll);
    else
        error('load_toronto_polygons:BadGeoJSON', ...
              'GeoJSON has no "features" or "coordinates" field.');
    end

    % Filter to plausible Toronto region (just defensive)
    mask = lonAll > -80 & lonAll < -78 & latAll > 43 & latAll < 45;
    lonAll = lonAll(mask);
    latAll = latAll(mask);

    if numel(lonAll) < 3
        warning('load_toronto_polygons:TooFewPoints', ...
                'Not enough points after filtering; returning empty.');
        polys = struct('lon', {}, 'lat', {});
        return;
    end

    % ---- Build a concave boundary polygon around the point cloud ----
    % shrinkFactor in (0,1]; smaller => more detailed / concave.
    shrinkFactor = 0.4;
    try
        k = boundary(lonAll, latAll, shrinkFactor);
    catch
        % Fallback to convex hull if boundary() not available
        k = convhull(lonAll, latAll);
    end

    polys(1).lon = lonAll(k);
    polys(1).lat = latAll(k);

    fprintf('load_toronto_polygons: built outline polygon with %d vertices.\n', numel(k));
end

% =====================================================================
% Recursive helper: walk any nesting of "coordinates" and extract [lon,lat]
% =====================================================================

function [lonAll, latAll] = collect_coords(coords, lonAll, latAll)
    if isnumeric(coords)
        % Whatever the original shape is, flatten to 2 columns
        c = coords;
        c = reshape(c, [], 2);   % each row is [lon, lat]
        lonAll = [lonAll; c(:,1)];
        latAll = [latAll; c(:,2)];
    elseif iscell(coords)
        for i = 1:numel(coords)
            [lonAll, latAll] = collect_coords(coords{i}, lonAll, latAll);
        end
    else
        % unknown type, ignore
    end
end
