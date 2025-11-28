%% run_dynamic_sim_map.m
% Phase 3: Dynamic Lloyd-style coverage with calls, battery, and bases,
%          visualised on top of a Toronto point cloud.
%
% Logic is the same as run_dynamic_sim.m (sim2). The only real change
% is the visualization: we convert [0,1]x[0,1] positions into
% longitude/latitude using a 16 km x 16 km bounding box that you chose
% with build_toronto_grid.m, and draw the density grid as a
% semi-transparent heatmap over Toronto.
%
% Assumptions:
%   - Phase 1 (run_lloyds_static.m) has already been run.
%   - The final base positions P (n x 2) are still in the workspace,
%     or stored as basePositions.
%   - build_toronto_grid.m has been run once, so that
%     data/toronto_area/toronto_grid_bbox.mat exists.
%   - static_density_map.m (sim1) and compute_coverage_centroids.m,
%     sample_call_location.m (sim2) are on the path.

clc; close all;
set(0, 'DefaultFigureWindowStyle', 'normal');   % open figures undocked

%% ---- Add paths to sim1 and sim2 helpers ----
thisFile   = mfilename('fullpath');
thisFolder = fileparts(thisFile);          % .../src/sim3
srcFolder  = fileparts(thisFolder);        % .../src
projRoot   = fileparts(srcFolder);         % .../Naloxone...

sim1Folder = fullfile(srcFolder, 'sim1');  % static Lloyd + density
sim2Folder = fullfile(srcFolder, 'sim2');  % dynamic helpers
addpath(sim1Folder, sim2Folder);

dataFolder = fullfile(projRoot, 'data', 'toronto_area');

%% ---- Check for base positions from Phase 1 ----
if exist('P','var')
    basePos = P;
elseif exist('basePositions','var')
    basePos = basePositions;
else
    error(['Base positions not found in workspace. ' ...
           'Run Phase 1 (run_lloyds_static.m) first so that P (n x 2) ' ...
           'is defined, then run this script.']);
end
nAgents = size(basePos, 1);

%% ---- Grid and density field (same as sim2, in [0,1]x[0,1]) ----
N  = 100;                % 100 x 100 grid (same as Phase 1)
xv = linspace(0, 1, N);
yv = linspace(0, 1, N);
[X, Y] = meshgrid(xv, yv);

D = static_density_map(X, Y);
if any(D(:) < 0)
    warning('Some densities are negative. Clipping to zero.');
    D = max(D, 0);
end
if all(D(:) == 0)
    error('Density map is all zeros. Cannot sample calls.');
end

%% ---- Load 16x16 km bounding box and Toronto point cloud ----
bboxFile = fullfile(dataFolder, 'toronto_grid_bbox.mat');
if ~isfile(bboxFile)
    error('run_dynamic_sim_map:NoBBox', ...
        ['Bounding box file not found at %s.\n' ...
         'Run build_toronto_grid.m first to choose the 16x16 km area.'], ...
         bboxFile);
end
S = load(bboxFile);
gridParams = S.gridParams;
lonMin = gridParams.lonMin;
lonMax = gridParams.lonMax;
latMin = gridParams.latMin;
latMax = gridParams.latMax;

% 1D lon/lat vectors corresponding to the 100x100 density grid
lonVec = lonMin + xv*(lonMax - lonMin);   % 1 x N
latVec = latMin + yv*(latMax - latMin);   % 1 x N

% Toronto point cloud from GeoJSON (already working)
shapes = load_toronto_geojson();
allLon = shapes.lon(:);
allLat = shapes.lat(:);

% Filter to plausible Toronto region
validMask = (allLon > -80 & allLon < -78) & ...
            (allLat >  43 & allLat <  45);
lonT = allLon(validMask);
latT = allLat(validMask);

if isempty(lonT)
    warning('run_dynamic_sim_map:EmptyCloud', ...
            'No filtered Toronto points found – background will be empty.');
end

%% ---- Physical / timing parameters (same as sim2) ----
domainSize_km = 16;          % 8 cells * 2 km each
cellSize_km   = domainSize_km / N;   % 0.16 km per 100x100 cell

droneSpeed_kmh   = 70;       % real drone speed (tweak as needed)
callRatePerHour  = 20;       % expected calls per hour (lambda)
serviceTime_min  = 5;        % time spent at call location
serviceTime_hr   = serviceTime_min / 60;

% Battery model
drainPercentPerHour    = 100;   % 100% per hour => 1hr continuous flight
chargePercentPerHour   = 150;   % >100% means full in <1h
batteryThresholdActive = 50;    % must be >= 50% to be dispatched / cover
reserveMarginPercent   = 5;     % extra percent to safely get home

% Simulation time and real-time scaling
simDuration_hr = 2.0;         % total simulated hours
dtSim_min      = 0.1;         % time step in minutes (0.1 min = 6 s sim time)
dtSim_hr       = dtSim_min / 60;
nSteps         = floor(simDuration_hr / dtSim_hr);

timeAccel = 100;              % sim runs (timeAccel)x faster than real time
pauseTime = dtSim_hr * 3600 / timeAccel;   % seconds to pause per loop

%% ---- Voronoi visualization toggle (same idea as sim2) ----
vorChoice = menu('Show Voronoi region map during simulation (normalised coords)?', ...
                 'Yes', ...
                 'No');
if vorChoice == 0
    error('User cancelled Voronoi visualization selection.');
end
showVoronoi = (vorChoice == 1);

%% ---- Drone state initialisation ----
STATE_IDLE_AT_BASE = 0;
STATE_TO_CALL      = 1;
STATE_AT_CALL      = 2;
STATE_COVER        = 3;
STATE_RETURNING    = 4;

pos      = basePos;                  % current positions in [0,1]^2
state    = STATE_IDLE_AT_BASE * ones(nAgents, 1);
battery  = 100 * ones(nAgents, 1);   % start fully charged
serviceEndTime = -inf * ones(nAgents, 1);
callAssigned   = nan(nAgents, 1);    % which call id each drone is on
coverageTarget = basePos;            % target for coverage / returning

% Call storage
calls = struct('pos', {}, 'startTime', {}, 'servedBy', {}, ...
               'done', {}, 'dispatchTime', {}, 'arrivalTime', {});

% For stats: actual travel times in minutes
arrivalTimesMin = [];

simTime_hr = 0;

%% ---- Figures ----
f1 = figure(1); clf(f1);
if showVoronoi
    f2 = figure(2); clf(f2);
    axVor = axes('Parent', f2);
end

%% ---- Main simulation loop ----
for step = 1:nSteps
    simTime_hr = simTime_hr + dtSim_hr;

    %% 1) Call arrival (Poisson) ---------------------------
    lambda = callRatePerHour;
    pCall  = lambda * dtSim_hr;
    if rand < pCall
        [cx, cy, ~, ~] = sample_call_location(X, Y, D);
        newCall.pos          = [cx, cy];
        newCall.startTime    = simTime_hr;
        newCall.servedBy     = NaN;
        newCall.done         = false;
        newCall.dispatchTime = NaN;
        newCall.arrivalTime  = NaN;
        calls(end+1) = newCall;  %#ok<SAGROW>

        callIdx = numel(calls);
        availableMask = (state == STATE_IDLE_AT_BASE | state == STATE_COVER) ...
                         & (battery >= batteryThresholdActive);
        if any(availableMask)
            idxAvail = find(availableMask);
            dists    = vecnorm(pos(idxAvail,:) - newCall.pos, 2, 2);
            [~, kMin] = min(dists);
            droneIdx  = idxAvail(kMin);

            state(droneIdx)        = STATE_TO_CALL;
            coverageTarget(droneIdx,:) = newCall.pos;
            callAssigned(droneIdx) = callIdx;
            calls(callIdx).servedBy    = droneIdx;
            calls(callIdx).dispatchTime = simTime_hr;

            dist_norm = norm(pos(droneIdx,:) - newCall.pos);
            dist_km   = dist_norm * domainSize_km;
            eta_hr    = dist_km / droneSpeed_kmh;
            eta_min   = eta_hr * 60;

            fprintf(['t = %.1f min: New call at (%.3f, %.3f) assigned to drone %d, ' ...
                     'ETA ~ %.1f min\n'], ...
                    simTime_hr*60, cx, cy, droneIdx, eta_min);
        else
            fprintf('t = %.1f min: New call at (%.3f, %.3f) but no drone available.\n', ...
                    simTime_hr*60, cx, cy);
        end
    end

    %% 2) Battery update -------------------------------
    for i = 1:nAgents
        if state(i) == STATE_IDLE_AT_BASE
            battery(i) = battery(i) + chargePercentPerHour * dtSim_hr;
        else
            battery(i) = battery(i) - drainPercentPerHour * dtSim_hr;
        end
        battery(i) = max(min(battery(i), 100), 0);
    end

    %% 3) Decide which drones must return on battery grounds -----
    for i = 1:nAgents
        if state(i) ~= STATE_IDLE_AT_BASE
            dHome_norm = hypot(pos(i,1) - basePos(i,1), pos(i,2) - basePos(i,2));
            dHome_km   = dHome_norm * domainSize_km;
            tReturn_hr = dHome_km / droneSpeed_kmh;

            batteryNeeded = drainPercentPerHour * tReturn_hr + reserveMarginPercent;

            if battery(i) <= batteryNeeded && state(i) ~= STATE_RETURNING
                fprintf('t = %.1f min: Drone %d low battery (%.1f%%). Returning to base.\n', ...
                        simTime_hr*60, i, battery(i));
                state(i) = STATE_RETURNING;
                coverageTarget(i,:) = basePos(i,:);
            end
        end
    end

    %% 4) Movement towards targets -------------------------
    maxStep_km = droneSpeed_kmh * dtSim_hr;

    for i = 1:nAgents
        switch state(i)
            case STATE_TO_CALL
                target = coverageTarget(i,:);
            case STATE_COVER
                target = coverageTarget(i,:);
            case STATE_RETURNING
                target = basePos(i,:);
            case STATE_AT_CALL
                target = pos(i,:);
            case STATE_IDLE_AT_BASE
                target = basePos(i,:);
            otherwise
                target = pos(i,:);
        end

        if state(i) == STATE_TO_CALL || state(i) == STATE_COVER || state(i) == STATE_RETURNING
            delta     = target - pos(i,:);
            dist_norm = norm(delta);
            dist_km   = dist_norm * domainSize_km;

            if dist_km <= maxStep_km || dist_norm == 0
                pos(i,:) = target;
            else
                stepFactor = maxStep_km / dist_km;
                pos(i,:) = pos(i,:) + stepFactor * delta;
            end
        end
    end

    %% 5) Arrival at call / service completion / returning ------
    for i = 1:nAgents
        switch state(i)
            case STATE_TO_CALL
                cIdx = callAssigned(i);
                if ~isnan(cIdx) && cIdx <= numel(calls)
                    cPos = calls(cIdx).pos;
                    if norm(pos(i,:) - cPos) * domainSize_km < 1e-3
                        state(i) = STATE_AT_CALL;
                        serviceEndTime(i) = simTime_hr + serviceTime_hr;

                        calls(cIdx).arrivalTime = simTime_hr;
                        if ~isnan(calls(cIdx).dispatchTime)
                            travelMin = (simTime_hr - calls(cIdx).dispatchTime) * 60;
                            arrivalTimesMin(end+1) = travelMin; %#ok<AGROW>
                        else
                            travelMin = NaN;
                        end

                        fprintf('t = %.1f min: Drone %d arrived at call %d. Travel time = %.1f min\n', ...
                                simTime_hr*60, i, cIdx, travelMin);
                    end
                end

            case STATE_AT_CALL
                if simTime_hr >= serviceEndTime(i)
                    cIdx = callAssigned(i);
                    if ~isnan(cIdx) && cIdx <= numel(calls)
                        calls(cIdx).done = true;
                    end
                    callAssigned(i) = NaN;

                    if battery(i) >= batteryThresholdActive + reserveMarginPercent
                        state(i) = STATE_COVER;
                        fprintf('t = %.1f min: Drone %d finished call, joining coverage.\n', ...
                                simTime_hr*60, i);
                    else
                        state(i) = STATE_RETURNING;
                        coverageTarget(i,:) = basePos(i,:);
                        fprintf('t = %.1f min: Drone %d finished call but low battery (%.1f%%). Returning.\n', ...
                                simTime_hr*60, i, battery(i));
                    end
                end

            case STATE_RETURNING
                if norm(pos(i,:) - basePos(i,:)) * domainSize_km < 1e-3
                    state(i) = STATE_IDLE_AT_BASE;
                    coverageTarget(i,:) = basePos(i,:);
                    fprintf('t = %.1f min: Drone %d arrived at base (battery %.1f%%).\n', ...
                            simTime_hr*60, i, battery(i));
                end
        end
    end

    %% 6) Lloyd coverage update for available drones -------------
    busyMask = (state == STATE_TO_CALL) | (state == STATE_AT_CALL) | (state == STATE_RETURNING);
    if any(busyMask)
        canCoverMask = (state == STATE_IDLE_AT_BASE | state == STATE_COVER) ...
                        & (battery >= batteryThresholdActive);
        coverageIdx = find(canCoverMask);

        if ~isempty(coverageIdx)
            P_cov = pos(coverageIdx,:);
            centroids_cov = compute_coverage_centroids(P_cov, X, Y, D);

            coverageTarget(coverageIdx,:) = centroids_cov;

            for k = 1:numel(coverageIdx)
                idx = coverageIdx(k);
                if state(idx) == STATE_IDLE_AT_BASE
                    state(idx) = STATE_COVER;
                end
            end
        end
    else
        for i = 1:nAgents
            if state(i) ~= STATE_IDLE_AT_BASE
                state(i) = STATE_RETURNING;
                coverageTarget(i,:) = basePos(i,:);
            end
        end
    end

    %% 7) Visualisation on Toronto map ---------------------------

    % Convert base and drone positions from [0,1]x[0,1] to lon/lat
    baseLon = lonMin + basePos(:,1) .* (lonMax - lonMin);
    baseLat = latMin + basePos(:,2) .* (latMax - latMin);

    droneLon = lonMin + pos(:,1) .* (lonMax - lonMin);
    droneLat = latMin + pos(:,2) .* (latMax - latMin);

    % Count active vs completed calls and compute avg arrival time
    if isempty(calls)
        nActiveCalls    = 0;
        nCompletedCalls = 0;
    else
        isDone     = [calls.done];
        isAssigned = ~isnan([calls.servedBy]);
        nActiveCalls    = sum(~isDone & isAssigned);
        nCompletedCalls = sum(isDone);
    end

    if isempty(arrivalTimesMin)
        avgArrivalStr = 'N/A';
    else
        avgArrivalStr = sprintf('%.1f', mean(arrivalTimesMin));
    end

    % Build list of active call positions in lon/lat
    activeCallLonLat = [];
    for c = 1:numel(calls)
        if ~calls(c).done && ~isnan(calls(c).servedBy)
            cPos = calls(c).pos;   % [x,y] in [0,1]
            cLon = lonMin + cPos(1) * (lonMax - lonMin);
            cLat = latMin + cPos(2) * (latMax - latMin);
            activeCallLonLat(end+1,:) = [cLon, cLat]; %#ok<SAGROW>
        end
    end

    % Main map figure
    set(0, 'CurrentFigure', f1);
    clf(f1);
    hold on;

    %% Background Toronto polygons
polys = load_toronto_polygons();
for k = 1:numel(polys)
    patch(polys(k).lon, polys(k).lat, [0.9 0.9 0.9], ...
          'EdgeColor', [0.4 0.4 0.4], 'LineWidth', 0.5);
end

    % Density grid as semi-transparent heatmap
    hImg = imagesc(lonVec, latVec, D);
    set(gca, 'YDir', 'normal');             % lat increases upwards
    set(hImg, 'AlphaData', 0.4);            % 0 = fully transparent, 1 = opaque

    axis equal tight;
    colorbar;
    xlabel('Longitude');
    ylabel('Latitude');

    % Active calls
    if ~isempty(activeCallLonLat)
        scatter(activeCallLonLat(:,1), activeCallLonLat(:,2), 80, 'm', 'p', 'filled', ...
                'DisplayName', 'Active calls');
    end

    % Bases
    scatter(baseLon, baseLat, 80, 'k', 's', 'filled', 'DisplayName', 'Bases');

    % Drones by state
    for i = 1:nAgents
        switch state(i)
            case STATE_IDLE_AT_BASE
                col = 'w'; mkr = 'o';
            case STATE_TO_CALL
                col = 'g'; mkr = 'o';
            case STATE_AT_CALL
                col = 'r'; mkr = 'o';
            case STATE_COVER
                col = 'c'; mkr = 'o';
            case STATE_RETURNING
                col = 'y'; mkr = 'o';
            otherwise
                col = 'w'; mkr = 'x';
        end

        scatter(droneLon(i), droneLat(i), 100, col, mkr, 'filled');
        txt = sprintf('%d (%.0f%%)', i, battery(i));
        text(droneLon(i), droneLat(i), [' ' txt], ...
             'Color','k', 'FontWeight','bold', 'FontSize',8);
    end

    title(sprintf(['Dynamic Lloyds on Toronto | t = %.1f min | ' ...
                   'active calls = %d | completed calls = %d | ' ...
                   'avg arrival = %s min'], ...
                  simTime_hr*60, nActiveCalls, nCompletedCalls, avgArrivalStr));

    hold off;
    drawnow;

    %% Optional: Voronoi region map (still in [0,1] coords) ------------
    if showVoronoi
        coverageMaskVor = (state == STATE_IDLE_AT_BASE | state == STATE_COVER) ...
                           & (battery >= batteryThresholdActive);
        coverageIdxVor = find(coverageMaskVor);

        cla(axVor);

        if isempty(coverageIdxVor)
            axis(axVor, 'off');
            text(axVor, 0.5, 0.5, 'No active coverage drones', ...
                 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'middle', ...
                 'FontWeight', 'bold');
        else
            xVecGrid = X(:);
            yVecGrid = Y(:);
            nCells   = numel(xVecGrid);
            owner    = zeros(nCells, 1);

            for k = 1:nCells
                dx = pos(coverageIdxVor,1) - xVecGrid(k);
                dy = pos(coverageIdxVor,2) - yVecGrid(k);
                distSq = dx.^2 + dy.^2;
                [~, localIdx] = min(distSq);
                owner(k) = coverageIdxVor(localIdx);
            end

            regionMapDyn = reshape(owner, size(X));

            imagesc(axVor, xv, yv, regionMapDyn);
            axis(axVor, 'xy', 'equal', 'tight');
            colorbar(axVor);
            hold(axVor, 'on');

            scatter(axVor, pos(coverageIdxVor,1), pos(coverageIdxVor,2), 80, 'k', 'filled');
            for idx = coverageIdxVor'
                text(axVor, pos(idx,1), pos(idx,2), sprintf(' %d', idx), ...
                     'Color','w', 'FontWeight','bold', 'FontSize',8);
            end

            title(axVor, sprintf('Dynamic Voronoi (coverage drones only) | t = %.1f min', ...
                                 simTime_hr*60));
            hold(axVor, 'off');
        end
    end

    pause(pauseTime);
end

fprintf('Dynamic Toronto-map simulation finished.\n');
