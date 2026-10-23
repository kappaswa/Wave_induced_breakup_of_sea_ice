clear; clc;

%% Files and geometry
initialFile   = 'PartFloating_0000.csv';
referenceFile = 'DE_connection_reference.mat';

deSize = [0.06, 0.06, 0.05];   % DE dimensions [m]
dp      = 0.005;                % particle spacing [m]
tol     = 1e-6;

% Include particles at dp/2 and 3dp/2 from an interface
faceBand = 1.5*dp + tol;

%% Read required columns
required = { ...
    'Pos.x [m]', ...
    'Pos.y [m]', ...
    'Pos.z [m]', ...
    'Idp', ...
    'DemBond'};

opts = detectImportOptions(initialFile, ...
    'Delimiter', ';', ...
    'NumHeaderLines', 3, ...
    'VariableNamingRule', 'preserve');

missing = setdiff(required, opts.VariableNames);

if ~isempty(missing)
    error('Missing CSV columns: %s', strjoin(missing, ', '));
end

opts.SelectedVariableNames = required;
T = readtable(initialFile, opts);

x       = double(T.('Pos.x [m]'));
y       = double(T.('Pos.y [m]'));
z       = double(T.('Pos.z [m]'));
idp     = uint32(T.Idp);
demBond = double(T.DemBond);

clear T

if any(~ismember(demBond, [0, 1, 2]))
    error('DemBond contains values other than 0, 1, or 2.');
end

%% Determine the DE grid from particle positions

xMin = min(x);
yMin = min(y);
zMin = min(z);

xMax = max(x);
yMax = max(y);
zMax = max(z);

origin = [xMin, yMin, zMin] - dp/2;

% Center-to-center span plus one particle spacing gives the full length.
numX = round((xMax - xMin + dp) / deSize(1));
numY = round((yMax - yMin + dp) / deSize(2));
numZ = round((zMax - zMin + dp) / deSize(3));

if any([numX, numY, numZ] < 1)
    error('Invalid DE-grid dimensions.');
end

numTotalDE = numX*numY*numZ;

if numTotalDE > double(intmax('uint16'))
    error('Too many DEs for uint16 storage.');
end

fprintf('Grid: %d x %d x %d DEs\n', ...
    numX, numY, numZ);

fprintf('Sheet size: %.6f x %.6f x %.6f m\n', ...
    numX*deSize(1), ...
    numY*deSize(2), ...
    numZ*deSize(3));

%% Keep only initially bonded particles
isBonded = demBond == 1;

if ~any(isBonded)
    error('No DemBond = 1 particles were found.');
end

xB   = x(isBonded);
yB   = y(isBonded);
zB   = z(isBonded);
idpB = idp(isBonded);

clear x y z idp demBond

%% Find the nearest x-grid boundary

xFaceIndex = round((xB - origin(1)) / deSize(1));

xFacePosition = ...
    origin(1) + xFaceIndex*deSize(1);

distX = abs(xB - xFacePosition);

isInternalXFace = ...
    xFaceIndex >= 1 & ...
    xFaceIndex <= numX - 1;

nearX = ...
    isInternalXFace & ...
    distX <= faceBand;

%% Find the nearest y-grid boundary

yFaceIndex = round((yB - origin(2)) / deSize(2));

yFacePosition = ...
    origin(2) + yFaceIndex*deSize(2);

distY = abs(yB - yFacePosition);

isInternalYFace = ...
    yFaceIndex >= 1 & ...
    yFaceIndex <= numY - 1;

nearY = ...
    isInternalYFace & ...
    distY <= faceBand;

%% Classify particles by interface
isCorner = nearX & nearY;

xOnly = nearX & ~nearY;
yOnly = nearY & ~nearX;

notNearInterface = ~(nearX | nearY);

if any(notNearInterface)

    badRows = find(notNearInterface);

    UnassignedParticles = table( ...
        idpB(badRows), ...
        xB(badRows), ...
        yB(badRows), ...
        zB(badRows), ...
        xFaceIndex(badRows), ...
        yFaceIndex(badRows), ...
        distX(badRows), ...
        distY(badRows), ...
        'VariableNames', { ...
            'Idp', ...
            'X', ...
            'Y', ...
            'Z', ...
            'NearestXFaceIndex', ...
            'NearestYFaceIndex', ...
            'DistanceToXFace', ...
            'DistanceToYFace'});

    fprintf('\nFirst unassigned particles:\n');

    disp(UnassignedParticles( ...
        1:min(20, height(UnassignedParticles)), :));

    writetable(UnassignedParticles, ...
        'unassigned_bond_particles.csv');

    error([ ...
        '%d DemBond = 1 particles are not near an internal ' ...
        'x- or y-interface. See unassigned_bond_particles.csv.'], ...
        nnz(notNearInterface));
end

fprintf('\nInterface classification\n');
fprintf('X-interface particles: %d\n', nnz(xOnly));
fprintf('Y-interface particles: %d\n', nnz(yOnly));
fprintf('Corner particles: %d\n', nnz(isCorner));

%% Orientation codes
orientation = zeros(size(idpB), 'uint8');

orientation(xOnly) = uint8(1);
orientation(yOnly) = uint8(2);

% 0 = corner
% 1 = x-interface
% 2 = y-interface

%% Initialize DE-pair arrays
de1 = zeros(size(idpB), 'uint16');
de2 = zeros(size(idpB), 'uint16');

%% Determine DE pairs for x-interface particles
%
% The interface is identified first.
%
% x-interface i separates:
%   left column  = i
%   right column = i + 1

if any(xOnly)

    leftColumn  = xFaceIndex(xOnly);
    rightColumn = leftColumn + 1;

    % Determine only the y-row and z-layer needed for this interface.
    rowY = floor( ...
        (yB(xOnly) - origin(2) + tol) / deSize(2)) + 1;

    if numZ == 1
        layerZ = ones(nnz(xOnly), 1);
    else
        layerZ = floor( ...
            (zB(xOnly) - origin(3) + tol) / deSize(3)) + 1;
    end

    rowY   = min(max(rowY,   1), numY);
    layerZ = min(max(layerZ, 1), numZ);

    leftDE = ...
        leftColumn ...
        + (rowY - 1)*numX ...
        + (layerZ - 1)*numX*numY;

    rightDE = ...
        rightColumn ...
        + (rowY - 1)*numX ...
        + (layerZ - 1)*numX*numY;

    de1(xOnly) = uint16(leftDE);
    de2(xOnly) = uint16(rightDE);
end

%% Determine DE pairs for y-interface particles
%
% y-interface j separates:
%   lower row = j
%   upper row = j + 1

if any(yOnly)

    lowerRow = yFaceIndex(yOnly);
    upperRow = lowerRow + 1;

    % Determine only the x-column and z-layer needed for this interface.
    columnX = floor( ...
        (xB(yOnly) - origin(1) + tol) / deSize(1)) + 1;

    if numZ == 1
        layerZ = ones(nnz(yOnly), 1);
    else
        layerZ = floor( ...
            (zB(yOnly) - origin(3) + tol) / deSize(3)) + 1;
    end

    columnX = min(max(columnX, 1), numX);
    layerZ  = min(max(layerZ,  1), numZ);

    lowerDE = ...
        columnX ...
        + (lowerRow - 1)*numX ...
        + (layerZ - 1)*numX*numY;

    upperDE = ...
        columnX ...
        + (upperRow - 1)*numX ...
        + (layerZ - 1)*numX*numY;

    de1(yOnly) = uint16(lowerDE);
    de2(yOnly) = uint16(upperDE);
end

%% Group particles belonging to the same interface
nonCorner = xOnly | yOnly;

keys = [ ...
    double(de1(nonCorner)), ...
    double(de2(nonCorner)), ...
    double(orientation(nonCorner))];

[uniqueKeys, ~, keyIndex] = unique(keys, 'rows');

numInterfaces = size(uniqueKeys, 1);

if numInterfaces > double(intmax('uint16'))
    error('Too many interfaces for uint16 storage.');
end

interfaceID = zeros(size(idpB), 'uint16');
interfaceID(nonCorner) = uint16(keyIndex);

numReferenceParticles = uint32( ...
    accumarray(keyIndex, 1, [numInterfaces, 1]));

%% Calculate containing-DE information for output compatibility
%
% These values are calculated after interface assignment.
% They are not used to identify the interface.

deXB = floor( ...
    (xB - origin(1) + tol) / deSize(1)) + 1;

deYB = floor( ...
    (yB - origin(2) + tol) / deSize(2)) + 1;

if numZ == 1
    deZB = ones(size(zB));
else
    deZB = floor( ...
        (zB - origin(3) + tol) / deSize(3)) + 1;
end

deXB = min(max(deXB, 1), numX);
deYB = min(max(deYB, 1), numY);
deZB = min(max(deZB, 1), numZ);

deIDB = deXB ...
      + (deYB - 1)*numX ...
      + (deZB - 1)*numX*numY;

%% Particle-level reference table
BondParticles = table( ...
    idpB, ...
    single(xB), ...
    single(yB), ...
    single(zB), ...
    uint16(deIDB), ...
    uint16(deXB), ...
    uint16(deYB), ...
    uint16(deZB), ...
    logical(isCorner), ...
    orientation, ...
    interfaceID, ...
    de1, ...
    de2, ...
    'VariableNames', { ...
        'Idp', ...
        'X0', ...
        'Y0', ...
        'Z0', ...
        'DE_ID', ...
        'DE_X', ...
        'DE_Y', ...
        'DE_Z', ...
        'IsCorner', ...
        'Orientation', ...
        'InterfaceID', ...
        'DE1', ...
        'DE2'});

%% Interface-level reference table
Interfaces = table( ...
    uint16((1:numInterfaces)'), ...
    uint16(uniqueKeys(:,1)), ...
    uint16(uniqueKeys(:,2)), ...
    uint8(uniqueKeys(:,3)), ...
    numReferenceParticles, ...
    'VariableNames', { ...
        'InterfaceID', ...
        'DE1', ...
        'DE2', ...
        'Orientation', ...
        'NumReferenceParticles'});

%% Check the expected number of interfaces
expectedXInterfaces = ...
    (numX - 1)*numY*numZ;

expectedYInterfaces = ...
    numX*(numY - 1)*numZ;

expectedInterfaces = ...
    expectedXInterfaces + expectedYInterfaces;

numFoundX = nnz(Interfaces.Orientation == 1);
numFoundY = nnz(Interfaces.Orientation == 2);

if height(Interfaces) ~= expectedInterfaces || ...
        numFoundX ~= expectedXInterfaces || ...
        numFoundY ~= expectedYInterfaces

    error([ ...
        'Interface grid is incomplete.\n' ...
        'Total: found %d, expected %d.\n' ...
        'X:     found %d, expected %d.\n' ...
        'Y:     found %d, expected %d.'], ...
        height(Interfaces), expectedInterfaces, ...
        numFoundX, expectedXInterfaces, ...
        numFoundY, expectedYInterfaces);
end

%% Save reference
Reference.InitialFile = initialFile;
Reference.ElementSize = deSize;
Reference.ParticleSpacing = dp;
Reference.FaceBand = faceBand;
Reference.Origin = origin;
Reference.NumDE = [numX, numY, numZ];

Reference.OrientationCodes = struct( ...
    'Corner', uint8(0), ...
    'X', uint8(1), ...
    'Y', uint8(2));

Reference.BondParticles = BondParticles;
Reference.Interfaces = Interfaces;

save(referenceFile, 'Reference', '-v7.3');

%% Summary
fprintf('\nReference file created: %s\n', referenceFile);
fprintf('Bonded particles: %d\n', numel(idpB));
fprintf('Corner particles: %d\n', nnz(isCorner));
fprintf('X-interface particles: %d\n', nnz(xOnly));
fprintf('Y-interface particles: %d\n', nnz(yOnly));
fprintf('Interfaces created: %d\n', height(Interfaces));
fprintf('All expected interfaces found: YES\n');

%% Top-view plot
fig = figure('Color', 'w');

scatter(xB, yB, 3, 'r', 'filled', ...
    'MarkerFaceAlpha', 0.85, ...
    'MarkerEdgeAlpha', 0);

hold on

sheetLx = numX*deSize(1);
sheetLy = numY*deSize(2);

rectangle( ...
    'Position', [origin(1), origin(2), sheetLx, sheetLy], ...
    'EdgeColor', 'k', ...
    'LineWidth', 2.5);

axis equal

xlim([origin(1), origin(1) + sheetLx]);
ylim([origin(2), origin(2) + sheetLy]);

xlabel('x [m]');
ylabel('y [m]');

title(sprintf( ...
    'Initial reference: %d x %d DEs', ...
    numX, numY), ...
    'FontWeight', 'bold');

set(gca, ...
    'FontSize', 14, ...
    'LineWidth', 1.2, ...
    'Box', 'on');

% exportgraphics(fig, ...
%     'reference_top_view.png', ...
%     'Resolution', 300);