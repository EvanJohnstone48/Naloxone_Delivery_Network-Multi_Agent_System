% run_lloyds_custom.m
clearvars; close all; clc;

% ==========================================
% PARAMETERS
% ==========================================
densityFile = '';        % '' to use density_map() without filename
N_AGENTS = 20;           % Number of agents
MAX_ITER = 100;           % Number of timesteps
FPS = 10;                % Animation speed

% Arena/Algorithm Parameters
ARENA_SIZE = 16;          % Coordinate bound (e.g., -3 to 3)
R_COM = 20.0;             % Radius of communication 
R_OBS = 6.0;             % Radius of observation 

% ==========================================
% 1. SETUP DENSITY MAP
% ==========================================
if ~isempty(densityFile) && exist(densityFile,'file')
    rho = readmatrix(densityFile);
else
    % Define the grid for sampling the density function
    M = 100; K = 100; % Resolution
    [xg, yg] = meshgrid(linspace(-ARENA_SIZE, ARENA_SIZE, K), ...
                        linspace(-ARENA_SIZE, ARENA_SIZE, M));
    t0 = 0;
    rho = density_map(xg, yg, t0);
end

% Normalize Density
rho = max(rho, 0);
if sum(rho(:)) > 0
    rho = rho / sum(rho(:));
else
    error('Density map is zero everywhere.');
end

% ==========================================
% 2. INITIALIZATION
% ==========================================
% P: (Time x Agents x 2) - Stores history of positions
P = nan(MAX_ITER, N_AGENTS, 2);

% G: (Agents x Agents x Time) - Stores history of adjacency
G = nan(N_AGENTS, N_AGENTS, MAX_ITER-1);

% Random Initial Positions in [-ARENA_SIZE, ARENA_SIZE]
init_pos = ARENA_SIZE * (2 * rand(N_AGENTS, 2) - 1);
P(1,:,:) = init_pos;

% Time vector
t = linspace(0, 10, MAX_ITER)'; 

% ==========================================
% 3. SIMULATION LOOP (Lloyd's Algorithm)
% ==========================================
fprintf('Starting Simulation (%d agents, %d steps)...\n', N_AGENTS, MAX_ITER);

for i = 1:MAX_ITER-1
    % Get current positions (N x 2)
    p_current = squeeze(P(i,:,:));
    
    % --- A. Calculate Adjacency (Connectivity) ---
    % This uses your src/lloyds_adjacency_matrix.m function correctly
    adj_matrix = lloyds_adjacency_matrix(p_current, R_COM);
    G(:,:,i) = adj_matrix;
    
    % --- B. Compute Centroids (Distributed) ---
    % 1. Find connected components (observation groups)
    graph_comp = conncomp(graph(adj_matrix, 'omitselfloops'));
    
    p_next = p_current; % Default: stay put if no move calculated
    
    for j = 1:max(graph_comp)
        % Indices of agents in this group
        group_idx = find(graph_comp == j);
        n_group = length(group_idx);
        p_group = p_current(group_idx, :);
        
        % 2. Get observed data points (x, y, density)
        % Uses your src/observation_set.m
        obs_data = observation_set(p_group, xg, yg, rho, R_OBS);
        
        if isempty(obs_data)
            continue; 
        end
        
        % 3. Voronoi Partition (K-Means)
        % Cluster observed points to nearest agent in the group
        if size(obs_data, 1) >= n_group
             [idx_cluster, ~] = kmeans(obs_data(:,1:2), n_group, ...
                'Start', p_group, 'EmptyAction', 'singleton', 'MaxIter', 10);
            
            % 4. Calculate Weighted Centroid for each agent
            for k = 1:n_group
                % Points belonging to k-th agent in this group
                cell_points = obs_data(idx_cluster == k, :);
                
                if ~isempty(cell_points)
                    mass = sum(cell_points(:,3)); % Sum of density
                    if mass > 0
                        cx = sum(cell_points(:,1) .* cell_points(:,3)) / mass;
                        cy = sum(cell_points(:,2) .* cell_points(:,3)) / mass;
                        
                        % Update position for next step
                        agent_global_id = group_idx(k);
                        p_next(agent_global_id, :) = [cx, cy];
                    end
                end
            end
        end
    end
    
    % Store next positions
    P(i+1,:,:) = p_next;
end

% ==========================================
% 4. PLOTTING & ANIMATION
% ==========================================
fprintf('Simulation Complete. Generating Plots...\n');

% Final Positions Plot
figure('Name','Final positions','NumberTitle','off');
imagesc(xg(1,:), yg(:,1), flipud(rho)); 
colormap(parula); set(gca, 'YDir', 'normal'); hold on;
axis square;
title('Final Node Positions over Density');

final_pos = squeeze(P(end,:,:));
scatter(final_pos(:,1), final_pos(:,2), 50, 'w', 'filled', 'MarkerEdgeColor', 'k');
xlim([-ARENA_SIZE ARENA_SIZE]);
ylim([-ARENA_SIZE ARENA_SIZE]);

% Run Animation if available
if exist('animate_network.m', 'file')
    try
        % animate_network.m is a script that expects t, P, G in workspace
        run('animate_network.m');
    catch ME
        warning(ME.identifier, 'Animation failed: %s', ME.message);
    end
else
    warning('animate_network.m not found. Skipping animation.');
end