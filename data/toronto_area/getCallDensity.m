% getCallDensity.m
% Estimate overdose-call "density" for a single grid square.
%
% Uses:
%   Call-weighted score:   D_c  = sum_i f_i * C_i
%   Area-adjusted density: rho_c = sum_i f_i * (C_i / A_i)
%
% where:
%   f_i = fraction of this square covered by neighbourhood i (0..1)
%   C_i = total calls in neighbourhood i
%   A_i = area of neighbourhood i in km^2


%% ---------------- Calls lookup: name -> calls ---------------------------

neighbourhoodNames = {
    'Downtown Yonge East'
    'Moss Park'
    'Kensington-Chinatown'
    'Church-Wellesley'
    'Yonge-Bay Corridor'
    'Harbourfront-CityPlace'
    'St Lawrence-East Bayfront-The Islands'
    'South Parkdale'
    'Fort York-Liberty Village'
    'West Queen West'
    'Kingsview Village-The Westway'
    'North St.James Town'
    'Wellington Place'
    'Trinity-Bellwoods'
    'Oakdale-Beverley Heights'
    'Annex'
    'East End-Danforth'
    'West Hill'
    'University'
    'Weston'
    'Junction-Wallace Emerson'
    'Oakridge'
    'Regent Park'
    'Tam O''Shanter-Sullivan'
    'New Toronto'
    'Dufferin Grove'
    'Bay-Cloverhill'
    'South Riverdale'
    'Islington'
    'Birchcliffe-Cliffside'
    'Forest Hill South'
    'Pelmo Park-Humberlea'
    'Golfdale-Cedarbrae-Woburn'
    'Caledonia-Fairbank'
    'Wexford/Maryvale'
    'High Park-Swansea'
    'Stonegate-Queensway'
    'Junction Area'
    'Corso Italia-Davenport'
    'Mount Dennis'
    'Weston-Pelham Park'
    'Palmerston-Little Italy'
    'Etobicoke City Centre'
    'South Eglinton-Davisville'
    'Henry Farm'
    'Agincourt South-Malvern West'
    'Eringate-Centennial-West Deane'
    'Ionview'
    'Keelesdale-Eglinton West'
    'Dorset Park'
    'Yonge-Eglinton'
    'Victoria Village'
    'Rexdale-Kipling'
    'Humber Summit'
    'Wychwood'
    'Runnymede-Bloor West Village'
    'Taylor-Massey'
    'Edenbridge-Humber Valley'
    'Oakwood Village'
    'Mimico-Queensway'
    'Little Portugal'
    'Rosedale-Moore Park'
    'Roncesvalles'
    'Woodbine Corridor'
    'Eglinton East'
    'Greenwood-Coxwell'
    'Briar Hill-Belgravia'
    'West Humber-Clairville'
    'Kennedy Park'
    'Clairlea-Birchmount'
    'High Park North'
    'Newtonbrook West'
    'Cabbagetown-South St.James Town'
    'North Riverdale'
    'The Beaches'
    'Rockcliffe-Smythe'
    'Playter Estates-Danforth'
    'Humewood-Cedarvale'
    'Woburn North'
    'York University Heights'
    'Danforth'
    'Glenfield-Jane Heights'
    'Downsview'
    'Yorkdale-Glen Park'
    'Black Creek'
    'Cliffcrest'
    'Willowridge-Martingrove-Richview'
    'Humber Bay Shores'
    'Elms-Old Rexdale'
    'Dovercourt Village'
    'Brookhaven-Amesbury'
    'Thistletown-Beaumond Heights'
    'Bendale-Glen Andrew'
    'Mount Olive-Silverstone-Jamestown'
};

neighbourhoodCalls = [
    445  % Downtown Yonge East
    424  % Moss Park
    156  % Kensington-Chinatown
    127  % Church-Wellesley
    121  % Yonge-Bay Corridor
    114  % Harbourfront-CityPlace
    102  % St Lawrence-East Bayfront-The Islands
    84   % South Parkdale
    80   % Fort York-Liberty Village
    71   % West Queen West
    57   % Kingsview Village-The Westway
    50   % North St.James Town
    47   % Wellington Place
    42   % Trinity-Bellwoods
    39   % Oakdale-Beverley Heights
    37   % Annex
    34   % East End-Danforth
    33   % West Hill
    32   % University
    32   % Weston
    31   % Junction-Wallace Emerson
    31   % Oakridge
    30   % Regent Park
    29   % Tam O''Shanter-Sullivan
    28   % New Toronto
    27   % Dufferin Grove
    24   % Bay-Cloverhill
    23   % South Riverdale
    23   % Islington
    9    % Birchcliffe-Cliffside
    9    % Forest Hill South
    9    % Pelmo Park-Humberlea
    8    % Golfdale-Cedarbrae-Woburn
    8    % Caledonia-Fairbank
    8    % Wexford/Maryvale
    8    % High Park-Swansea
    8    % Stonegate-Queensway
    8    % Junction Area
    8    % Corso Italia-Davenport
    8    % Mount Dennis
    7    % Weston-Pelham Park
    7    % Palmerston-Little Italy
    7    % Etobicoke City Centre
    7    % South Eglinton-Davisville
    7    % Henry Farm
    7    % Agincourt South-Malvern West
    7    % Eringate-Centennial-West Deane
    7    % Ionview
    7    % Keelesdale-Eglinton West
    7    % Dorset Park
    6    % Yonge-Eglinton
    6    % Victoria Village
    6    % Rexdale-Kipling
    6    % Humber Summit
    6    % Wychwood
    5    % Runnymede-Bloor West Village
    5    % Taylor-Massey
    5    % Edenbridge-Humber Valley
    21   % Oakwood Village
    20   % Mimico-Queensway
    19   % Little Portugal
    18   % Rosedale-Moore Park
    17   % Roncesvalles
    17   % Woodbine Corridor
    16   % Eglinton East
    16   % Greenwood-Coxwell
    15   % Briar Hill-Belgravia
    14   % West Humber-Clairville
    14   % Kennedy Park
    14   % Clairlea-Birchmount
    14   % High Park North
    13   % Newtonbrook West
    13   % Cabbagetown-South St.James Town
    13   % North Riverdale
    12   % The Beaches
    12   % Rockcliffe-Smythe
    12   % Playter Estates-Danforth
    11   % Humewood-Cedarvale
    11   % Woburn North
    11   % York University Heights
    11   % Danforth
    11   % Glenfield-Jane Heights
    11   % Downsview
    10   % Yorkdale-Glen Park
    10   % Black Creek
    5    % Cliffcrest
    5    % Willowridge-Martingrove-Richview
    5    % Humber Bay Shores
    5    % Elms-Old Rexdale
    5    % Dovercourt Village
    5    % Brookhaven-Amesbury
    5    % Thistletown-Beaumond Heights
    5    % Bendale-Glen Andrew
    5    % Mount Olive-Silverstone-Jamestown
];

callKeys = lower(strtrim(neighbourhoodNames));
callMap  = containers.Map(callKeys, neighbourhoodCalls);

%% ---------------- Area lookup: name -> km^2 -----------------------------

scriptFolder = fileparts(mfilename('fullpath'));
csvPathAreas = fullfile(scriptFolder,'toronto_areas.csv');
areaTable    = readtable(csvPathAreas);

areaKeys = lower(strtrim(areaTable.name));
areaVals = areaTable.area_km2;
areaMap  = containers.Map(areaKeys, areaVals);

%% ---------------- Number lookup: number -> name -------------------------
% CSV should have columns: number, name

csvPathNums = fullfile(scriptFolder,'toronto_neighbourhood_numbers.csv');
numTable    = readtable(csvPathNums);

idKeys = numTable.number;                    % numeric ids
idVals = lower(strtrim(numTable.name));      % neighbourhood names
idToName = containers.Map(idKeys, idVals);   % number -> lowercase name

%% ---------------- User input -------------------------------------------

n = input('How many neighbourhoods overlap this square? ');

callScore = 0;
densScore = 0;
missingArea = {};

for i = 1:n
    fprintf('\nNeighbourhood %d:\n', i);
    id  = input('  Neighbourhood number: ');
    pct = input('  Coverage of square (0–100): ');
    f   = pct / 100;
    
    if ~isKey(idToName, id)
        error('Neighbourhood number %d not found in toronto_neighbourhood_numbers.csv.', id);
    end
    
    nameKey = idToName(id);  % lowercase, trimmed
    fprintf('    → %s\n', nameKey);  % optional: echo resolved name
    
    % --- calls ---
    if isKey(callMap, nameKey)
        C = callMap(nameKey);
    else
        fprintf('    (No call data for "%s" → using C = 0)\n', nameKey);
        C = 0;
    end
    
    callScore = callScore + f * C;
    
    % --- area (if available) ---
    if isKey(areaMap, nameKey)
        A = areaMap(nameKey);
        densScore = densScore + f * (C / A);
    else
        missingArea{end+1} = nameKey; %#ok<AGROW>
    end
end

%% ---------------- Output -----------------------------------------------

fprintf('\n====================================\n');
fprintf('Call-weighted score:       %.2f\n', callScore);
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
%fprintf('Water/white regions: don''t enter them; that part of the square\n');
%fprintf('implicitly contributes 0 calls and 0 density.\n');