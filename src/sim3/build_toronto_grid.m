function build_toronto_grid()
% build_toronto_grid  Interactively place a 16 km x 16 km grid on Toronto.
%
% Uses neighbourhood polygons from toronto.geojson as background.
% Saves bounding box as data/toronto_area/toronto_grid_bbox.mat

    clc; close all;

    % Load polygons
    polys = load_toronto_polygons();
    if isempty(polys)
        error('build_toronto_grid:NoPolys', ...
              'No usable polygons found in toronto.geojson.');
    end

    % Gather all coords for axis limits
    allLon = [];
    allLat = [];
    for k = 1:numel(polys)
        allLon = [allLon; polys(k).lon(:)];
        allLat = [allLat; polys(k).lat(:)];
    end

    % Plot neighbourhood polygons
    f = figure('Name','Toronto neighbourhood polygons', ...
               'NumberTitle','off');
    axes('Parent', f);
    hold on;

    for k = 1:numel(polys)
        % Light grey fill, darker edges
        patch(polys(k).lon, polys(k).lat, [0.9 0.9 0.9], ...
              'EdgeColor', [0.4 0.4 0.4], 'LineWidth', 0.5);
    end

    axis equal;
    xlim([min(allLon) max(allLon)]);
    ylim([min(allLat) max(allLat)]);
    grid on;
    xlabel('Longitude');
    ylabel('Latitude');
    title({'Toronto neighbourhood polygons', ...
           'Click BOTTOM-LEFT then TOP-RIGHT corners of your 16 km x 16 km grid'});

    % --- User clicks bottom-left and top-right corners ---
    disp('Click bottom-left corner of the 16 km x 16 km grid...');
    [xb, yb] = ginput(1);
    plot(xb, yb, 'rx', 'MarkerSize', 10, 'LineWidth', 2);

    disp('Click top-right corner of the 16 km x 16 km grid...');
    [xt, yt] = ginput(1);
    plot(xt, yt, 'rx', 'MarkerSize', 10, 'LineWidth', 2);

    lonMin = min(xb, xt);
    lonMax = max(xb, xt);
    latMin = min(yb, yt);
    latMax = max(yb, yt);

    % Approximate physical size of the box (km) just for sanity check
    R_earth = 6371;  % km
    midLatRad   = mean([latMin latMax]) * pi/180;
    kmPerDegLat = pi * R_earth / 180;
    kmPerDegLon = kmPerDegLat * cos(midLatRad);

    dx_km = (lonMax - lonMin) * kmPerDegLon;
    dy_km = (latMax - latMin) * kmPerDegLat;

    fprintf('Selected box approx %.2f km (lon) by %.2f km (lat)\n', dx_km, dy_km);

    % Draw the rectangle for visual confirmation
    rectangle('Position', [lonMin, latMin, lonMax - lonMin, latMax - latMin], ...
              'EdgeColor', 'r', 'LineWidth', 1.5);
    drawnow;

    % Ask if we should save
    choice = questdlg( ...
        sprintf('Save this box? (%.1f km x %.1f km)', dx_km, dy_km), ...
        'Confirm grid', ...
        'Yes', 'No', 'Yes');

    if ~strcmp(choice, 'Yes')
        disp('User cancelled; nothing saved.');
        return;
    end

    % Prepare output struct
    gridParams.lonMin = lonMin;
    gridParams.lonMax = lonMax;
    gridParams.latMin = latMin;
    gridParams.latMax = latMax;
    gridParams.dx_km  = dx_km;
    gridParams.dy_km  = dy_km;

    % Save to data/toronto_area/toronto_grid_bbox.mat
    thisFile   = mfilename('fullpath');
    sim3Folder = fileparts(thisFile);
    srcFolder  = fileparts(sim3Folder);
    projRoot   = fileparts(srcFolder);
    dataFolder = fullfile(projRoot, 'data', 'toronto_area');
    outFile    = fullfile(dataFolder, 'toronto_grid_bbox.mat');

    save(outFile, 'gridParams');
    fprintf('Saved grid bounding box to %s\n', outFile);
end
