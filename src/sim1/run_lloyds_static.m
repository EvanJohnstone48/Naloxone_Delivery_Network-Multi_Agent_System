clear; clc; close all; % Clear workspace, command window, and close all figures
set(0, 'DefaultFigureWindowStyle', 'normal');

%% Make the grid for high resolution density
N = 100;                      % 100x100 density grid
xv = linspace(0, 1, N); % this just makes 100 x points between the coordinate system 0 and 1, returns a vector
yv = linspace(0, 1, N); % Same as x
[X, Y] = meshgrid(xv, yv); % Creates the 2D grid from these vectors

clear static_density_map   %clear this unless we dont need to reinitialize density map
D = static_density_map(X, Y);     % Input our X and Y matrices and get the density map

% Basic sanity check
if any(D(:) < 0)
    warning('Some densities are negative. Clipping to zero.');
    D = max(D, 0);
end

%% Choose number of drones
defaultN = 4;
userN = input(sprintf('Enter number of drones [default = %d]: ', defaultN));
if isempty(userN)
    nAgents = defaultN;
else
    nAgents = userN;
end
fprintf('Using %d drones.\n', nAgents)

%% Select Lloyd update mode (pure vs relaxed)
lloydMode = menu('Select Lloyd update mode:', ...
                 'Pure Lloyd (alpha = 1)', ...
                 'Relaxed Lloyd (choose alpha in (0,1))');
if lloydMode == 0
    error('User cancelled Lloyd mode selection.');
end

switch lloydMode
    case 1  % Pure Lloyd
        alpha = 1.0;
        fprintf('Using pure Lloyd: alpha = 1.0 (teleport to centroids each iter).\n');
    case 2  % Relaxed Lloyd
        defaultAlpha = 0.3;
        alphaInput = input(sprintf('Enter alpha in (0,1) [default = %.2f]: ', defaultAlpha));
        if isempty(alphaInput)
            alpha = defaultAlpha;
        else
            alpha = alphaInput;
        end
        if alpha <= 0 || alpha > 1
            error('alpha must be in (0,1]. You entered %.3f', alpha);
        end
        fprintf('Using relaxed Lloyd: alpha = %.3f (fractional move towards centroids).\n', alpha);
end

%% Toggle for centroid markers/arrows
centroidChoice = menu('Show centroid markers/arrows on density plot?', ...
                      'Yes', ...
                      'No');
if centroidChoice == 0
    error('User cancelled centroid visualization selection.');
end
showCentroids = (centroidChoice == 1);


%% Choose initialization mode
initChoice = menu('How to initialise drone positions?', ...
                  'Random in domain [0,1]x[0,1]', ...
                  'Click positions on density map');
if initChoice == 0
    error('Initialization choice cancelled by user.');
end

% Initialize positions P as nAgents x 2 matrix
P = zeros(nAgents, 2);

switch initChoice
    case 1  % Random
        rng('shuffle'); %sets a random seed for rand function
        P(:, 1) = rand(nAgents, 1);  % x in [0,1], this function conviently matches our coordinate domain
        P(:, 2) = rand(nAgents, 1);  % y in [0,1]
        fprintf('Initialised drone positions randomly.\n');
    case 2  % Click
        figure('Name','Click initial positions','NumberTitle','off', 'DefaultFigureWindowStyle', 'normal');
        imagesc(xv, yv, D); 
        axis xy equal tight;
        colorbar;
        title(sprintf('Click %d initial drone positions', nAgents)); % Update title dynamically
        xlabel('x');
        ylabel('y');
        hold on;
        
        fprintf('Please click %d points on the figure for initial positions.\n', nAgents);
        
        % --- CHANGED: Use a loop to plot points as you click ---
        xClick = zeros(nAgents, 1); % Pre-allocate arrays
        yClick = zeros(nAgents, 1);
        
        for i = 1:nAgents
            % Get exactly 1 click
            [xi, yi] = ginput(1); 
            
            % Store the coordinate
            xClick(i) = xi;
            yClick(i) = yi;
            
            % Plot a red filled dot immediately
            plot(xi, yi, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
            
            % Force MATLAB to draw it NOW (otherwise it might wait until the loop ends)
            drawnow; 
        end
        % -------------------------------------------------------
        
        hold off;
        pause(0.5); % Optional: Short pause so you can see your final clicks
        close;
        
        P(:, 1) = xClick(:);
        P(:, 2) = yClick(:);
        fprintf('Initialised drone positions from user clicks.\n');
end

%% Lloyd iteration parameters
maxIter = 100;         % maximum number of iterations
tol     = 1e-3;       % stop if max movement smaller than this
pauseTime = 0.1;      % pause between frames for visualization
J_history = nan(maxIter, 1);  % to store Lloyd cost per iteration, for our graph

fprintf('\nStarting Lloyd iterations...\n');
fprintf('alpha = %.2f, maxIter = %d, tol = %.1e\n\n', alpha, maxIter, tol);

% ----- Create figures once so we don't keep stealing focus -----
f1 = figure(1);        % main visualization (density + Voronoi)
clf(f1); %Clear contents

f2 = figure(2);        % cost plot figure
clf(f2);
axCost = axes('Parent', f2);   % axes for J vs iteration, added it here so it isnt constantly updating

% ---------------- MODIFIED SECTION 1: LOAD IMAGE ----------------
% Added try/catch to prevent crash if file is missing
% Added flipud because 'axis xy' puts (0,0) at bottom, but images put (0,0) at top
try
    backgroundMap = imread('../../data/Density_Map_Images/cropped_neighbourhood_map.jpg');
    backgroundMap = flipud(backgroundMap); 
    fprintf('Background map loaded successfully.\n');
catch
    warning('Background map not found. Proceeding with white background.');
    backgroundMap = [];
end
% ----------------------------------------------------------------

%% Main Lloyd loop
for iter = 1:maxIter
    % One Lloyd step
    [P_new, centroids, regionMap, maxMove, J] = lloyd_step(P, X, Y, D, alpha);
    J_history(iter) = J;

    % Print some info to the command window
    fprintf('Iteration %2d: max move = %.4f\n', iter, maxMove);
    for j = 1:nAgents
        fprintf('  Drone %d centroid = (%.3f, %.3f), pos = (%.3f, %.3f)\n', ...
                j, centroids(j,1), centroids(j,2), P_new(j,1), P_new(j,2));
    end
    fprintf('\n');

    % Plot current state
    figure(1); clf;

    % Left: density heatmap
    subplot(1, 2, 1); %Divides into 2 subplots, makes 1st subplot active

    % ---------------- MODIFIED SECTION 2: PLOT LAYERS ----------------
    % 1. Draw the Background Map first (Bottom Layer)
    if ~isempty(backgroundMap)
        image([0 1], [0 1], backgroundMap); 
        hold on;
    end
    
    % 2. Draw the Density Map ON TOP (Middle Layer)
    hDensity = imagesc(xv, yv, D);
    
    % 3. Set Transparency (Only if we have a background map)
    if ~isempty(backgroundMap)
        set(hDensity, 'AlphaData', 0.6); % 0.6 = 60% opaque
    end
    
    % Formatting
    axis xy equal tight;
    colorbar;
    title(sprintf('Density map with drones (iter %d)', iter));
    xlabel('x');
    ylabel('y');
    hold on;
    % -----------------------------------------------------------------

    %Plot the drones stored in P_new on top of imagem 80 is marker size,
    %'w' is marker face coloour, 'MarkerEdgeColor','k' is border colour
    scatter(P_new(:,1), P_new(:,2), 80, 'w', 'filled', 'MarkerEdgeColor','k');

    % Label drones
    for j = 1:nAgents
        %text anchored at p_new location
        text(P_new(j,1), P_new(j,2), sprintf(' %d', j), ...
            'Color','w', 'FontWeight','bold', 'FontSize',10);
    end

    % Optional: centroid markers + arrows showing where drones "want" to go
    if showCentroids
        % Centroid markers (red crosses)
        scatter(centroids(:,1), centroids(:,2), 60, 'r', 'x', 'LineWidth', 1.5);

        % Arrows/segments from drone to its centroid
        for j = 1:nAgents
            plot([P_new(j,1), centroids(j,1)], ...
                 [P_new(j,2), centroids(j,2)], ...
                 'w--', 'LineWidth', 1.0);
        end
    end

    hold off;

    % Right: Voronoi regions on grid
    subplot(1, 2, 2);
    imagesc(xv, yv, regionMap); %Regionmap already has all the region data from lloyds step
    axis xy equal tight;
    colorbar;
    title('Voronoi regions (grid based)');
    xlabel('x');
    ylabel('y');
    hold on;
    scatter(P_new(:,1), P_new(:,2), 80, 'k', 'filled');
    for j = 1:nAgents
        text(P_new(j,1), P_new(j,2), sprintf(' %d', j), ...
            'Color','w', 'FontWeight','bold', 'FontSize',10);
    end
    hold off;

    drawnow; %updates the figure windows :)

    %%Figure 2: Lloyd cost J vs iteration (update without stealing focus)
    cla(axCost);   % clear only the J-plot axes
    iterVec = 1:iter;
    plot(axCost, iterVec, J_history(1:iter), '-o');
    xlabel(axCost, 'Iteration');
    ylabel(axCost, 'Cost J');
    title(axCost, 'Lloyd cost vs iteration');
    grid(axCost, 'on');

    drawnow;
    pause(pauseTime);

    % Check stopping condition
    P = P_new;  % update for next iteration
    if maxMove < tol
        fprintf('Converged at iteration %d (max move < %.1e).\n', iter, tol);
        break;
    end
end

fprintf('Lloyd static simulation finished.\n');
%TOO EASYYYYY
