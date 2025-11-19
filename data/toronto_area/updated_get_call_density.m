% getCallDensity.m
% Estimate overdose-call "density" for a single grid square.

%% ---------------- Script folder -------------------------
scriptFolder = fileparts(mfilename('fullpath'));

%% ---------------- Load calls + names + IDs -------------------------
csvPathCalls = fullfile(scriptFolder, "neighbourhood_name_num_call.csv");
callTable    = readtable(csvPathCalls);

callNums  = callTable.("NeighbourhoodNumber");
callNames = lower(strtrim(callTable.("NeighbourhoodName")));
callVals  = callTable.("NumberOfCalls");

% Maps
idToCallName = containers.Map(callNums, callNames);
callMap      = containers.Map(callNames, callVals);

%% ---------------- Load area data (UPDATED) -------------------------
% Use the new fixed_neighborhood_areas.csv and look up by NAME
csvPathAreas = fullfile(scriptFolder, "toronto_neighbourhood_areas.csv");
areaTable    = readtable(csvPathAreas);

% Assuming columns: neighbourhood_number, neighbourhood_name, area_km2
% We ONLY care about the name + area; numbers in this file are not trusted.
areaKeys = lower(strtrim(areaTable.neighbourhood_name));
areaVals = areaTable.area_km2;

areaMap = containers.Map(areaKeys, areaVals);

%% ---------------- User input ----------------------------
n = input('How many neighbourhoods overlap this square? ');

callScore   = 0;
densScore   = 0;
missingArea = {};
totalPct    = 0;

for i = 1:n
    fprintf('\nNeighbourhood %d:\n', i);
    id  = input('  Neighbourhood number: ');
    pct = input('  Coverage of square (0–100): ');
    totalPct = totalPct + pct;

    f = pct / 100;

    %% ---- ID → name (from call table) ----
    try
        nameKey = idToCallName(id);
        fprintf('    → %s\n', nameKey);
    catch
        fprintf('    (Neighbourhood number %d not found → using name="" and C=0)\n', id);
        nameKey = "";
    end

    %% ---- Calls lookup by name ----
    if nameKey == ""
        C = 0;
    else
        try
            C = callMap(nameKey);
            fprintf('    Number of calls: %d\n', C);
        catch
            fprintf('    (No call data for "%s" → using Number of calls = 0)\n', nameKey);
            C = 0;
        end
    end

    %% ---- Area lookup by *name* only (numbers in area file ignored) ----
    if nameKey ~= ""
        try
            A = areaMap(nameKey);
            densScore = densScore + f * (C / A);
        catch
            missingArea{end+1} = nameKey; %#ok<AGROW>
        end
    end
end

%% ---------------- Validate percentages sum to 100 --------------------
if abs(totalPct - 100) > 1e-6
    fprintf('\nERROR: Total coverage = %.1f%% (must add up to 100%%).\n', totalPct);
    fprintf('Please restart and enter correct percentages.\n\n');
    return;
end

%% ---------------- Output -------------------------------------------
fprintf('\n====================================\n');

if densScore > 0
    fprintf('Area-adjusted density:     %.4f calls/km^2\n', densScore);
else
    fprintf('Area-adjusted density:     not available (no matching areas)\n');
end

if ~isempty(missingArea)
    fprintf('\nMissing area data for:\n');
    disp(unique(missingArea)');
end

fprintf('====================================\n');
