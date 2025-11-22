function d = static_density_map(x, y)

% this returns a smooth high-resolution density map
% Inputs:
%   x, y, n by n meshgrids of coordinates in [0,1] x [0,1]
%
% Output is d a n by n matrix of nonzero densities

    persistent baseX baseY baseD initialised %if we dont clear this will stay initialized

    if isempty(initialised) || isempty(baseD)
        % Ask user which density map to load
        choice = menu('Select density dataset:', ...
                      'Real-world map (Filled_in_heat_map)', ...
                      'Testing map (Filled_in_heat_map_testing)');
        if choice == 0
            error('density_map:NoSelection', ...
                  'No dataset selected. User canceled the menu.');
        end

        % This file is in .../src/sim2
        simFolder = fileparts(mfilename('fullpath'));  % .../src/sim2
        srcFolder = fileparts(simFolder);              % .../src
        projRoot  = fileparts(srcFolder);              % .../Naloxone...

        dataFolder = fullfile(projRoot, 'data', 'toronto_area');

        %Lets the user choose between csv files
        switch choice
            case 1
                csvName = 'Filled_in_heat_map.csv';
            case 2
                csvName = 'Filled_in_heat_map_testing.csv';
            otherwise
                error('Unexpected menu choice.');
        end

        csvPath = fullfile(dataFolder, csvName);
        fprintf('Loading density data from: %s\n', csvPath);

        baseD = readmatrix(csvPath);

        if isempty(baseD)
            error('density_map:EmptyData', ...
                  'File %s appears to be empty or unreadable.', csvPath);
        end

        % Ensure it is numeric and nonnegative
        baseD = double(baseD);
        %Finds every negative element and sets it to 0
        baseD(baseD < 0) = 0;

        [nRows, nCols] = size(baseD);
        fprintf('Loaded base density of size %d x %d.\n', nRows, nCols);

        % Build base coordinate grid in [0,1] x [0,1]
        [baseX, baseY] = meshgrid(linspace(0, 1, nCols), ...
                                  linspace(0, 1, nRows));

        initialised = true;
    end

    % Interpolate the 8x8 grid onto the 100x100 grid
    d = interp2(baseX, baseY, baseD, x, y, 'spline');

    % interp2 can sometimes return something call NaN?? Just make it 0
    d(isnan(d)) = 0;

    %normalisation so all densities are between 0 and 1
    maxVal = max(d(:));
    if maxVal > 0
    d = d / maxVal;
    end
end
