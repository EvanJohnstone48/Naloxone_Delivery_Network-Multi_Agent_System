%% run_dynamic_sim.m
% Phase 2: Dynamic Lloyd-style coverage with calls, battery, and bases.
%
% States:
%   0 = IDLE_AT_BASE
%   1 = TO_CALL
%   2 = AT_CALL
%   3 = COVER   (off base, doing Lloyd coverage)
%   4 = RETURNING
%

clc; close all;
set(0, 'DefaultFigureWindowStyle', 'normal');   % open figures undocked
thisFile   = mfilename('fullpath');
thisFolder = fileparts(thisFile);        % .../src/sim2
srcFolder  = fileparts(thisFolder);      % .../src
sim1Folder = fullfile(srcFolder, 'sim1');% .../src/sim1

addpath(sim1Folder);
fprintf('Added static sim folder to path: %s\n', sim1Folder);

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

%% ---- Grid and density field ----
N  = 100;                % 100 x 100 grid (same as Phase 1)
xv = linspace(0, 1, N);
yv = linspace(0, 1, N);
[X, Y] = meshgrid(xv, yv);


% Use the same density loader as Phase 1.
% If your function is named density_map, switch to that.
D = static_density_map(X, Y);
if any(D(:) < 0)
    warning('Some densities are negative. Clipping to zero.');
    D = max(D, 0);
end
if all(D(:) == 0)
    error('Density map is all zeros. Cannot sample calls.');
end

%% ---- Physical / timing parameters ----
domainSize_km = 16;          % 8 cells * 2 km each
cellSize_km   = domainSize_km / N;   % 0.16 km per 100x100 cell

droneSpeed_kmh = 70;         % real drone speed (tweak as needed)
callRatePerHour = 20;         % expected calls per hour (lambda)
serviceTime_min = 5;         % time spent at call location
serviceTime_hr  = serviceTime_min / 60;

% Battery mode
drainPercentPerHour   = 100;   % 100% per hour => 1hr continuous flight from full to empty
chargePercentPerHour  = 150;   % >100% means full charge in <1h
batteryThresholdActive = 50;   % must be >= 50% to be dispatched or cover
reserveMarginPercent   = 5;    % extra percent above minimum needed to return

% Simulation time and real-time scaling
simDuration_hr = 2.0;         % total simulated hours
dtSim_min      = 0.1;         % time step in minutes (0.5 min = 30 s sim time)
dtSim_hr       = dtSim_min / 60;
nSteps         = floor(simDuration_hr / dtSim_hr);

timeAccel = 1000;                % sim runs (n)x faster than real time
pauseTime = dtSim_hr * 3600 / timeAccel;   % seconds to pause per step

%% ---- Voronoi visualization toggle ----
vorChoice = menu('Show Voronoi region map during simulation?', ...
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

pos      = basePos;                  % current positions
state    = STATE_IDLE_AT_BASE * ones(nAgents, 1);
battery  = 100 * ones(nAgents, 1);   % start fully charged
serviceEndTime = -inf * ones(nAgents, 1);   % when AT_CALL finishes
callAssigned   = nan(nAgents, 1);    % which call id each drone is on
coverageTarget = basePos;            % current target for coverage / returning

% Call storage, dynamic list
calls = struct('pos', {}, 'startTime', {}, 'servedBy', {}, ...
               'done', {}, 'dispatchTime', {}, 'arrivalTime', {});

% For stats: actual travel times (arrival - dispatch) in minutes
arrivalTimesMin = [];

simTime_hr = 0;

%% ---- Figures ----
f1 = figure(1);
clf(f1);

% Optional Voronoi figure
if showVoronoi
    f2 = figure(2);
    clf(f2);
    axVor = axes('Parent', f2);
end

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

%% ---- Main simulation loop ----
for step = 1:nSteps
    simTime_hr = simTime_hr + dtSim_hr;

    %% 1) Call arrival (Poisson) ----------------------------------------
    lambda = callRatePerHour;
    pCall = lambda * dtSim_hr;    % prob of at least one call in this small interval
    if rand < pCall
        % Sample a call location from density map
        [cx, cy, ~, ~] = sample_call_location(X, Y, D);
        newCall.pos          = [cx, cy];
        newCall.startTime    = simTime_hr;
        newCall.servedBy     = NaN;
        newCall.done         = false;
        newCall.dispatchTime = NaN;
        newCall.arrivalTime  = NaN;
        calls(end+1) = newCall;  %#ok<SAGROW>

        % Assign nearest available drone
        callIdx = numel(calls);
        availableMask = (state == STATE_IDLE_AT_BASE | state == STATE_COVER) ...
                         & (battery >= batteryThresholdActive);
        if any(availableMask)
            idxAvail = find(availableMask);
            dists = vecnorm(pos(idxAvail,:) - newCall.pos, 2, 2);
            [~, kMin] = min(dists);
            droneIdx = idxAvail(kMin);

            % Dispatch that drone
            state(droneIdx)        = STATE_TO_CALL;
            coverageTarget(droneIdx,:) = newCall.pos;
            callAssigned(droneIdx) = callIdx;
            calls(callIdx).servedBy = droneIdx;
            calls(callIdx).dispatchTime = simTime_hr;

            % Compute and print theoretical ETA (real world minutes)
            dist_norm = norm(pos(droneIdx,:) - newCall.pos);   % in [0,1]
            dist_km   = dist_norm * domainSize_km;
            eta_hr    = dist_km / droneSpeed_kmh;
            eta_min   = eta_hr * 60;

            fprintf(['t = %.1f min: New call at (%.3f, %.3f) assigned to drone %d, ' ...
                 'ETA ~ %.1f min\n'], ...
                simTime_hr*60, cx, cy, droneIdx, eta_min);
        else
            fprintf('t = %.1f min: New call at (%.3f, %.3f) but no drone available (call unserved).\n', ...
                simTime_hr*60, cx, cy);
        end
    end

    %% 2) Battery update (drain / charge) -------------------------------
    for i = 1:nAgents
        if state(i) == STATE_IDLE_AT_BASE
            % Recharge at base
            battery(i) = battery(i) + chargePercentPerHour * dtSim_hr;
        else
            % Drain while off base / doing work
            battery(i) = battery(i) - drainPercentPerHour * dtSim_hr;
        end
        battery(i) = max(min(battery(i), 100), 0);  % clamp to [0,100]
    end

    %% 3) Decide which drones must return on battery grounds -----------
    for i = 1:nAgents
        if state(i) ~= STATE_IDLE_AT_BASE
            dHome_norm = hypot(pos(i,1) - basePos(i,1), pos(i,2) - basePos(i,2));
            dHome_km   = dHome_norm * domainSize_km;
            tReturn_hr = dHome_km / droneSpeed_kmh;

            batteryNeeded = drainPercentPerHour * tReturn_hr + reserveMarginPercent;

            if battery(i) <= batteryNeeded && state(i) ~= STATE_RETURNING
                % Flip this drone to RETURNING; it will head back now
                fprintf('t = %.1f min: Drone %d low battery (%.1f%%). Returning to base.\n', ...
                        simTime_hr*60, i, battery(i));
                state(i) = STATE_RETURNING;
                coverageTarget(i,:) = basePos(i,:);
            end
        end
    end

    %% 4) Update movement towards targets -------------------------------
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
                % Stay parked at call location; no movement
                target = pos(i,:);
            case STATE_IDLE_AT_BASE
                % Stay at base
                target = basePos(i,:);
            otherwise
                target = pos(i,:);
        end

        if state(i) == STATE_TO_CALL || state(i) == STATE_COVER || state(i) == STATE_RETURNING
            delta      = target - pos(i,:);
            dist_norm  = norm(delta);
            dist_km    = dist_norm * domainSize_km;

            if dist_km <= maxStep_km || dist_norm == 0
                % Arrive at target in this step
                pos(i,:) = target;
            else
                % Move along the direction with limited step
                stepFactor = (maxStep_km / dist_km);
                pos(i,:) = pos(i,:) + stepFactor * delta;
            end
        end
    end

    %% 5) Handle arrival at call and service completion ----------------
    for i = 1:nAgents
        switch state(i)
            case STATE_TO_CALL
            % If drone has reached its call location, flip to AT_CALL
            cIdx = callAssigned(i);
            if ~isnan(cIdx) && cIdx <= numel(calls)
                cPos = calls(cIdx).pos;
                if norm(pos(i,:) - cPos) * domainSize_km < 1e-3  % effectively there
                    state(i) = STATE_AT_CALL;
                    serviceEndTime(i) = simTime_hr + serviceTime_hr;

                    % Record arrival time and actual travel time
                    calls(cIdx).arrivalTime = simTime_hr;
                    if ~isnan(calls(cIdx).dispatchTime)
                        travelMin = (simTime_hr - calls(cIdx).dispatchTime) * 60;
                        arrivalTimesMin(end+1) = travelMin;  %#ok<AGROW>
                    else
                        travelMin = NaN;
                    end

                    fprintf('t = %.1f min: Drone %d arrived at call %d. Travel time = %.1f min\n', ...
                            simTime_hr*60, i, cIdx, travelMin);
                end
            end


            case STATE_AT_CALL
                % When service time is done, either join coverage or return
                if simTime_hr >= serviceEndTime(i)
                    cIdx = callAssigned(i);
                    if ~isnan(cIdx) && cIdx <= numel(calls)
                        calls(cIdx).done = true;
                    end
                    callAssigned(i) = NaN;

                    % Decide behaviour based on battery
                    if battery(i) >= batteryThresholdActive + reserveMarginPercent
                        state(i) = STATE_COVER;   % join coverage
                        fprintf('t = %.1f min: Drone %d finished call, joining coverage.\n', ...
                                simTime_hr*60, i);
                    else
                        state(i)        = STATE_RETURNING;
                        coverageTarget(i,:) = basePos(i,:);
                        fprintf('t = %.1f min: Drone %d finished call but low battery (%.1f%%). Returning.\n', ...
                                simTime_hr*60, i, battery(i));
                    end
                end

            case STATE_RETURNING
                % If reached base, flip to IDLE_AT_BASE
                if norm(pos(i,:) - basePos(i,:)) * domainSize_km < 1e-3
                    state(i) = STATE_IDLE_AT_BASE;
                    coverageTarget(i,:) = basePos(i,:);
                    fprintf('t = %.1f min: Drone %d arrived at base (battery %.1f%%).\n', ...
                            simTime_hr*60, i, battery(i));
                end
        end
    end

    %% 6) Lloyd coverage update for available drones -------------------
    % Only do coverage if at least one drone is off-line (on call or returning)
    busyMask = (state == STATE_TO_CALL) | (state == STATE_AT_CALL) | (state == STATE_RETURNING);
    if any(busyMask)
        % Drones that can cover: either at base or already in COVER, with enough battery
        canCoverMask = (state == STATE_IDLE_AT_BASE | state == STATE_COVER) ...
                        & (battery >= batteryThresholdActive);
        coverageIdx = find(canCoverMask);

        if ~isempty(coverageIdx)
            P_cov = pos(coverageIdx,:);
            centroids_cov = compute_coverage_centroids(P_cov, X, Y, D);

            % Update coverage targets
            coverageTarget(coverageIdx,:) = centroids_cov;

            % Any at-base drones that started covering now change state
            for k = 1:numel(coverageIdx)
                idx = coverageIdx(k);
                if state(idx) == STATE_IDLE_AT_BASE
                    state(idx) = STATE_COVER;
                end
            end
        end
    else
        % No busy drones: drift all drones back to base
        for i = 1:nAgents
            if state(i) ~= STATE_IDLE_AT_BASE
                state(i) = STATE_RETURNING;
                coverageTarget(i,:) = basePos(i,:);
            end
        end
    end

    %% 7) Visualisation ------------------------------------------------
    set(0, 'CurrentFigure', f1);
    clf(f1);
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % <--- CHANGED START: Layered Plotting with Transparency --->
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 1. Draw the Background Map first (Bottom Layer)
    if ~isempty(backgroundMap)
        % 'image' draws from top-left by default, but we flipped data.
        % This maps the image to the [0,1] x [0,1] domain.
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
    xlabel('x');
    ylabel('y');
    hold on;
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % <--- CHANGED END --->
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Plot call positions
    activeCallPos = [];
    for c = 1:numel(calls)
        if ~calls(c).done && ~isnan(calls(c).servedBy)
            activeCallPos(end+1,:) = calls(c).pos; %#ok<SAGROW>
        end
    end
    if ~isempty(activeCallPos)
        scatter(activeCallPos(:,1), activeCallPos(:,2), 80, 'm', 'p', 'filled', ...
                'DisplayName', 'Active calls');
    end
    
    % Count active vs completed calls for the title
    if isempty(calls)
        nActiveCalls    = 0;
        nCompletedCalls = 0;
    else
        isDone     = [calls.done];
        isAssigned = ~isnan([calls.servedBy]);
        % "Active" = not done AND assigned to a drone
        nActiveCalls    = sum(~isDone & isAssigned);
        nCompletedCalls = sum(isDone);
    end

    % Average arrival time (minutes) from arrivalTimesMin
    if isempty(arrivalTimesMin)
        avgArrivalStr = 'N/A';
    else
        avgArrivalStr = sprintf('%.1f', mean(arrivalTimesMin));
    end

    title(sprintf(['Dynamic Lloyds | t = %.1f min | active calls = %d | ' ...
               'completed calls = %d | avg arrival = %s min'], ...
               simTime_hr*60, nActiveCalls, nCompletedCalls, avgArrivalStr));
    
    % Plot base positions
    scatter(basePos(:,1), basePos(:,2), 80, 'k', 's', 'filled', ...
            'DisplayName', 'Bases');

    % Plot drones with different markers by state
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

        scatter(pos(i,1), pos(i,2), 100, col, mkr, 'filled');
        txt = sprintf('%d (%.0f%%)', i, battery(i));
        text(pos(i,1), pos(i,2), [' ' txt], ...
             'Color','k', 'FontWeight','bold', 'FontSize',8);
    end

    hold off;

    %% Optional: Voronoi region map in second figure -------------------
    if showVoronoi
        % Coverage-eligible drones: same idea as in Lloyd coverage step:
        %   - either at base or in COVER
        %   - with enough battery to be considered active
        coverageMaskVor = (state == STATE_IDLE_AT_BASE | state == STATE_COVER) ...
                           & (battery >= batteryThresholdActive);
        coverageIdxVor = find(coverageMaskVor);

        cla(axVor);

        if isempty(coverageIdxVor)
            % No active coverage drones: just show a message
            axis(axVor, 'off');
            text(axVor, 0.5, 0.5, 'No active coverage drones', ...
                 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'middle', ...
                 'FontWeight', 'bold');
        else
            % Compute grid-based Voronoi regions for COVERAGE drones only
            xVec = X(:);
            yVec = Y(:);
            nCells = numel(xVec);
            owner  = zeros(nCells, 1);

            for k = 1:nCells
                dx = pos(coverageIdxVor,1) - xVec(k);
                dy = pos(coverageIdxVor,2) - yVec(k);
                distSq = dx.^2 + dy.^2;
                [~, localIdx] = min(distSq);
                % Map back to GLOBAL drone index
                owner(k) = coverageIdxVor(localIdx);
            end

            regionMapDyn = reshape(owner, size(X));

            % Plot into the Voronoi axes without stealing focus
            imagesc(axVor, xv, yv, regionMapDyn);
            axis(axVor, 'xy', 'equal', 'tight');
            colorbar(axVor);
            hold(axVor, 'on');

            % Plot COVERAGE drones on top (the ones that own the regions)
            scatter(axVor, pos(coverageIdxVor,1), pos(coverageIdxVor,2), 80, 'k', 'filled');
            for idx = coverageIdxVor'
                text(axVor, pos(idx,1), pos(idx,2), sprintf(' %d', idx), ...
                     'Color','w', 'FontWeight','bold', 'FontSize',8);
            end

            title(axVor, sprintf('Dynamic Voronoi (coverage drones only) | t = %.1f min', simTime_hr*60));
            hold(axVor, 'off');
        end
    end

    drawnow;
    pause(pauseTime);
end

fprintf('Dynamic simulation finished.\n');
