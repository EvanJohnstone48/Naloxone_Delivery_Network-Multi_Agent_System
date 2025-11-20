
% plot_density_from_excel.m
% Read csv density grid, turn it into a matrix,
% then plot as 3D surface + heatmap.


choice = menu("Select Data Mode", "Testing", "RealWorld");
if choice == 1
    mode = "Testing";
elseif choice == 2
    mode = "RealWorld";
else
    error("No selection made.");
end

switch mode
    case "Testing"
        filename = "Filled_in_heat_map_testing.csv";
    case "RealWorld"
        filename = "Filled_in_heat_map.csv";
    otherwise
        error("Unknown mode: %s", mode);
end

fprintf("Mode: %s\nUsing file: %s\n", mode, filename);

%% --- READ DATA -------------------------------------------------------
M = readmatrix(filename);
disp("Loaded density matrix M (size):");
disp(size(M));


%%3D
[nRows, nCols] = size(M);
[x, y] = meshgrid(1:nCols, 1:nRows);

figure;
surf(x, y, M);
shading interp;          % smooth color between points
colormap hot;            % heat-style colormap
colorbar;
xlabel("Column index");
ylabel("Row index");
zlabel("Density");
title("3D Density Surface from Excel Grid");


% Plot the heatmap
%% --- HEATMAP: 0 = BLUE (water), >0 = smooth 50-color red palette -------

figure;

% Replace exact zeros with NaN so they map to MissingDataColor (blue)
Mplot = M;
Mplot(Mplot == 0) = NaN;

h = heatmap(Mplot);
xlabel("Column");
ylabel("Row");
title("Density Heatmap (0 = Blue, Smooth Red Scale)");

% Blue for water
h.MissingDataColor = [1 1 1];
h.MissingDataLabel = 'Zero';

% Base 5-tone red palette from your scale
baseRed = [
    1.00 0.97 0.96   % lightest
    0.99 0.91 0.89
    0.97 0.81 0.80
    0.94 0.66 0.66
    0.84 0.60 0.60   % darkest
];

% Interpolate to create a smooth ~50-color red colormap
numColors = 50;
xBase = linspace(0, 1, size(baseRed,1));
xInterp = linspace(0, 1, numColors);

smoothRed = interp1(xBase, baseRed, xInterp);

% Apply the 50-color red colormap
colormap(h, smoothRed);

% Color range only for NON-zero values
minNonZero = min(Mplot(:), [], 'omitnan');
maxVal     = max(Mplot(:), [], 'omitnan');
h.ColorLimits = [minNonZero maxVal];

h.ColorbarVisible = 'on';