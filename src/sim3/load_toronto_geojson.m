function shapes = load_toronto_geojson()
% load_toronto_geojson  Load Toronto coords from a GeoJSON file as a point cloud.
%
% Output:
%   shapes : struct with fields
%       shapes.lon  : column vector of longitudes
%       shapes.lat  : column vector of latitudes
%
% We don't try to reconstruct exact polygon rings. We just pull out
% all [lon,lat] pairs from every feature so we can plot Toronto.

    % Figure out data/toronto_area relative to this file
    thisFile   = mfilename('fullpath');
    thisFolder = fileparts(thisFile);          % .../src/sim3
    srcFolder  = fileparts(thisFolder);        % .../src
    projRoot   = fileparts(srcFolder);         % .../Naloxone...
    dataFolder = fullfile(projRoot, 'data', 'toronto_area');

    geojsonPath = fullfile(dataFolder, 'toronto.geojson');  % change name if needed
    if ~isfile(geojsonPath)
        error('load_toronto_geojson:FileNotFound', ...
              'GeoJSON file not found at: %s', geojsonPath);
    end

    fprintf('Loading GeoJSON from: %s\n', geojsonPath);
    txt = fileread(geojsonPath);
    gj  = jsondecode(txt);

    if ~isfield(gj, 'features')
        error('load_toronto_geojson:BadFormat', ...
              'GeoJSON file does not have a "features" field.');
    end

    features = gj.features;

    allLon = [];
    allLat = [];

    for i = 1:numel(features)
        if ~isfield(features(i),'geometry') || isempty(features(i).geometry)
            continue;
        end
        geom = features(i).geometry;
        if ~isfield(geom,'coordinates')
            continue;
        end

        % Recursively gather all numeric values from "coordinates"
        vals = gather_numbers(geom.coordinates);   % column vector

        if numel(vals) < 2
            continue;
        end

        % Interpret as [lon,lat,lon,lat,...]
        if mod(numel(vals), 2) ~= 0
            vals = vals(1:2*floor(numel(vals)/2));   % drop last odd element
        end

        lon_i = vals(1:2:end);
        lat_i = vals(2:2:end);

        allLon = [allLon; lon_i(:)];
        allLat = [allLat; lat_i(:)];
    end

    if isempty(allLon)
        warning('load_toronto_geojson:NoCoords', ...
                'No coordinate pairs found in GeoJSON. Plot will be empty.');
    end

    shapes = struct('lon', allLon, 'lat', allLat);

    fprintf('load_toronto_geojson: collected %d coordinate pairs in total.\n', ...
            numel(allLon));
end


function out = gather_numbers(x)
% gather_numbers  Recursively pull all numeric values out of x into a column vector.

    if isnumeric(x)
        out = x(:);
    elseif iscell(x)
        out = [];
        for k = 1:numel(x)
            out = [out; gather_numbers(x{k})]; %#ok<AGROW>
        end
    elseif isstruct(x)
        out = [];
        fn = fieldnames(x);
        for k = 1:numel(fn)
            out = [out; gather_numbers(x.(fn{k}))]; %#ok<AGROW>
        end
    else
        out = [];
    end
end

