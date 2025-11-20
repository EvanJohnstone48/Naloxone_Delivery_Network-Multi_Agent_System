function d = density_map(x, y, t)
%{
Inputs:
x, y: n by n meshgrid of x and y coordinates
t: time (unused in this static map, but required by interface)
Outputs:
d: n by n meshgrid of densities matching x,y resolution
%}

    % 1. Robust Relative Path (as discussed before)
    [currentFunctionPath, ~, ~] = fileparts(mfilename('fullpath'));
    csvFile = fullfile(currentFunctionPath, '..', 'data', 'toronto_area', 'Filled_in_heat_map.csv');

    if ~isfile(csvFile)
        error('File not found: %s', csvFile);
    end

    % 2. Read the raw density matrix
    rho_raw = readmatrix(csvFile);
    
    % 3. Create a coordinate system for the CSV data
    % We assume the CSV covers the full extent of the input x and y.
    % We create a temporary grid for the CSV so we can interpolate from it.
    
    % Get dimensions of the CSV
    [rows_raw, cols_raw] = size(rho_raw);
    
    % Create vectors spanning the min and max of the user's input area
    % (This maps the CSV to the physical bounds of x and y)
    x_vec = linspace(min(x(:)), max(x(:)), cols_raw);
    y_vec = linspace(min(y(:)), max(y(:)), rows_raw);
    
    [X_raw, Y_raw] = meshgrid(x_vec, y_vec);

    % 4. Interpolate 'rho_raw' onto the user's 'x, y' grid
    % 'linear' is standard; 'spline' or 'cubic' are smoother but slower.
    % We transpose rho_raw if necessary, but usually readmatrix matches meshgrid orientation.
    d = interp2(X_raw, Y_raw, rho_raw, x, y, 'linear', 0); 
    % Note: The '0' at the end sets any points outside the CSV range to 0.

    % 5. Normalize (ensure density constraints)
    % Ensure no negative values (interp2 can sometimes cause slight negatives)
    d(d < 0) = 0;
    
    % Normalize so sum is 1 (optional, depending on definition of density)
    if sum(d(:)) > 0
        d = d / sum(d(:));
    end
end