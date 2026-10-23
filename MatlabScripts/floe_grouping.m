clear;
clc;
close all;

%% Files
analysisFile = 'DE_remaining_bonds_0600.mat';
resultFile   = 'DE_floes_0600.mat';

%% Load analysis and reference
S = load(analysisFile, 'Analysis');
Analysis = S.Analysis;

Interfaces = Analysis.Interfaces;

R = load(Analysis.ReferenceFile, 'Reference');
Reference = R.Reference;

numX = double(Reference.NumDE(1));
numY = double(Reference.NumDE(2));
numZ = double(Reference.NumDE(3));

nDE = numX * numY * numZ;
numIfaces = height(Interfaces);

if numZ ~= 1
    error('This top-view script currently expects NumDE(3) = 1.');
end

%% Interface data
de1 = double(Interfaces.DE1);
de2 = double(Interfaces.DE2);

connectedInterface = logical(Interfaces.IsConnected);
brokenInterface    = ~connectedInterface;

%% Validate DE IDs
invalidID = ...
    de1 < 1 | de1 > nDE | ...
    de2 < 1 | de2 > nDE;

if any(invalidID)
    error('%d interfaces contain invalid DE IDs.', nnz(invalidID));
end

if any(de1 == de2)
    error('%d interfaces connect a DE to itself.', nnz(de1 == de2));
end

%% DE coordinates
DE_ID = (1:nDE)';

[DE_X, DE_Y, DE_Z] = ind2sub( ...
    [numX, numY, numZ], DE_ID);

DE_X = DE_X(:);
DE_Y = DE_Y(:);
DE_Z = DE_Z(:);

%% Check that interfaces connect neighboring DEs
dX = abs(DE_X(de1) - DE_X(de2));
dY = abs(DE_Y(de1) - DE_Y(de2));
dZ = abs(DE_Z(de1) - DE_Z(de2));

nonNeighborInterface = (dX + dY + dZ) ~= 1;

fprintf('\nINTERFACE CHECK\n');
fprintf('Total interfaces:              %d\n', numIfaces);
fprintf('Connected interfaces:          %d\n', nnz(connectedInterface));
fprintf('Broken interfaces:             %d\n', nnz(brokenInterface));
fprintf('Non-neighboring DE interfaces: %d\n', nnz(nonNeighborInterface));

if any(nonNeighborInterface)
    warning(['Some interfaces connect non-neighboring DEs. ', ...
             'The reference DE numbering should be checked.']);
end

%% Check for repeated DE pairs
pairs = sort([de1, de2], 2);
[uniquePairs, ~, pairGroup] = unique(pairs, 'rows');

rowsPerPair = accumarray(pairGroup, 1);
numRepeatedPairs = nnz(rowsPerPair > 1);

fprintf('Repeated DE pairs:             %d\n', numRepeatedPairs);

if numRepeatedPairs > 0
    warning('Some DE pairs appear more than once in Interfaces.');
end

%% Build graph from connected interfaces
s = de1(connectedInterface);
t = de2(connectedInterface);

G = graph(s, t, [], nDE);

%% Connected graph components are floes
floeID = conncomp(G);
floeID = floeID(:);

numFloes = max(floeID);

%% Classify broken interfaces
sameFloe = floeID(de1) == floeID(de2);

% Broken interface whose DEs remain connected through another path
brokenButBypassed = brokenInterface & sameFloe;

% Broken interface whose DEs belong to different components
actualFloeBoundary = brokenInterface & ~sameFloe;

fprintf('\nFLOE CONNECTIVITY\n');
fprintf('Broken but bypassed:           %d\n', nnz(brokenButBypassed));
fprintf('Actual floe-boundary segments: %d\n', nnz(actualFloeBoundary));
fprintf('Isolated DEs:                  %d\n', nnz(degree(G) == 0));
fprintf('Number of floes:               %d\n', numFloes);

%% Floe statistics
numDEPerFloe = accumarray( ...
    floeID, ...
    1, ...
    [numFloes, 1], ...
    @sum, ...
    0);

deArea = ...
    Reference.ElementSize(1) * ...
    Reference.ElementSize(2);

deVolume = ...
    Reference.ElementSize(1) * ...
    Reference.ElementSize(2) * ...
    Reference.ElementSize(3);

floeArea   = numDEPerFloe * deArea;
floeVolume = numDEPerFloe * deVolume;

FloeTable = table( ...
    (1:numFloes)', ...
    uint32(numDEPerFloe), ...
    floeArea, ...
    floeVolume, ...
    'VariableNames', { ...
        'FloeID', ...
        'NumDE', ...
        'Area_m2', ...
        'Volume_m3'});

DETable = table( ...
    uint32(DE_ID), ...
    uint16(DE_X), ...
    uint16(DE_Y), ...
    uint16(DE_Z), ...
    uint32(floeID), ...
    'VariableNames', { ...
        'DE_ID', ...
        'DE_X', ...
        'DE_Y', ...
        'DE_Z', ...
        'FloeID'});

fprintf('Largest floe:                 %d DEs\n', max(numDEPerFloe));
fprintf('Smallest floe:                %d DEs\n\n', min(numDEPerFloe));

%% Save results
Floes.ReferenceFile = Analysis.ReferenceFile;
Floes.AnalysisFile  = analysisFile;

Floes.Graph = G;

Floes.DETable   = DETable;
Floes.FloeTable = FloeTable;
Floes.FloeID    = floeID;

Floes.ConnectedInterface = connectedInterface;
Floes.BrokenInterface = brokenInterface;
Floes.BrokenButBypassed = brokenButBypassed;
Floes.ActualFloeBoundary = actualFloeBoundary;

Floes.NumFloes = numFloes;

save(resultFile, 'Floes', '-v7.3');

%% Plot colored floes in DE-index coordinates
gridFloe = reshape(floeID, [numX, numY])';

figure('Color', 'w');

imagesc(1:numX, 1:numY, gridFloe);

axis equal tight
set(gca, 'YDir', 'normal');

xlabel('DE index in x');
ylabel('DE index in y');

title(sprintf( ...
    'Floes from remaining interface graph: %d floes', ...
    numFloes), ...
    'FontWeight', 'bold');

colormap(turbo(max(numFloes, 2)));
clim([0.5, numFloes + 0.5]);

cb = colorbar;
cb.Label.String = 'Floe ID';

set(gca, ...
    'FontSize', 14, ...
    'LineWidth', 1.2, ...
    'Box', 'on');

drawnow

%% Physical DE coordinates
x0 = double(Reference.Origin(1));
y0 = double(Reference.Origin(2));

dxDE = double(Reference.ElementSize(1));
dyDE = double(Reference.ElementSize(2));

xCenter = x0 + (DE_X - 0.5) * dxDE;
yCenter = y0 + (DE_Y - 0.5) * dyDE;

sheetLx = numX * dxDE;
sheetLy = numY * dyDE;

%% Plot actual graph-derived floe boundaries
boundaryRows = find(actualFloeBoundary);

%% Compare top 5 mm broken bonds, fully broken interfaces,
%  and actual floe boundaries

particleStatus = Analysis.ParticleStatus;

%% Broken non-corner bond particles

nonCornerParticle = ~logical(particleStatus.IsCorner);

brokenParticle = ...
    double(particleStatus.DemBondFinal) == 0 & ...
    nonCornerParticle;

%% Select broken particles in the top 5 mm

zRef = double(particleStatus.Z0);

bandThickness = 5e-3;   % top 5 mm
zTol = 1e-10;

zMax = max(zRef(nonCornerParticle));

topBand = ...
    nonCornerParticle & ...
    zRef >= zMax - bandThickness - zTol;

topBrokenParticle = ...
    brokenParticle & topBand;

xTopBroken = ...
    double(particleStatus.X0(topBrokenParticle));

yTopBroken = ...
    double(particleStatus.Y0(topBrokenParticle));

fprintf('\nFIGURE 2 DATA\n');
fprintf('Top-band z range: %.6f to %.6f m\n', ...
    zMax - bandThickness, zMax);

fprintf('Broken particles in top 5 mm: %d\n', ...
    nnz(topBrokenParticle));

fprintf('Fully broken DE interfaces:   %d\n', ...
    nnz(brokenInterface));

fprintf('Actual floe boundaries:       %d\n', ...
    nnz(actualFloeBoundary));

%% Create Figure 2

figure('Color', 'w');
hold on

%% Grey: broken bond particles in the top 5 mm

hTopBroken = scatter( ...
    xTopBroken, ...
    yTopBroken, ...
    4, ...
    [0.0 0.0 0.55], ...
    'filled', ...
    'MarkerFaceAlpha', 0.15, ...
    'MarkerEdgeColor', 'none');

%% Red legend handle for fully broken interfaces

hBrokenInterface = plot( ...
    nan, ...
    nan, ...
    'r-', ...
    'LineWidth', 1.3);

%% Red: all fully broken DE interfaces

allBrokenRows = find(brokenInterface);

for q = 1:numel(allBrokenRows)

    k = allBrokenRows(q);

    a = de1(k);
    b = de2(k);

    xa = xCenter(a);
    ya = yCenter(a);

    xb = xCenter(b);
    yb = yCenter(b);

    xInterface = 0.5 * (xa + xb);
    yInterface = 0.5 * (ya + yb);

    % DEs differ in x:
    % interface is vertical in the top view
    if DE_X(a) ~= DE_X(b)

        plot( ...
            [xInterface, xInterface], ...
            [yInterface-dyDE/2, yInterface+dyDE/2], ...
            'r-', ...
            'LineWidth', 2.3, ...
            'HandleVisibility', 'off');

    % DEs differ in y:
    % interface is horizontal in the top view
    elseif DE_Y(a) ~= DE_Y(b)

        plot( ...
            [xInterface-dxDE/2, xInterface+dxDE/2], ...
            [yInterface, yInterface], ...
            'r-', ...
            'LineWidth', 2.3, ...
            'HandleVisibility', 'off');

    end
end

%% Black legend handle for actual floe boundaries

hFloeBoundary = plot( ...
    nan, ...
    nan, ...
    'k-', ...
    'LineWidth', 2.5);

%% Black: broken interfaces that separate different floes

floeBoundaryRows = find(actualFloeBoundary);

for q = 1:numel(floeBoundaryRows)

    k = floeBoundaryRows(q);

    a = de1(k);
    b = de2(k);

    xa = xCenter(a);
    ya = yCenter(a);

    xb = xCenter(b);
    yb = yCenter(b);

    xInterface = 0.5 * (xa + xb);
    yInterface = 0.5 * (ya + yb);

    % DEs differ in x:
    % interface is vertical in the top view
    if DE_X(a) ~= DE_X(b)

        plot( ...
            [xInterface, xInterface], ...
            [yInterface-dyDE/2, yInterface+dyDE/2], ...
            'k-', ...
            'LineWidth', 1.3);

    % DEs differ in y:
    % interface is horizontal in the top view
    elseif DE_Y(a) ~= DE_Y(b)

        plot( ...
            [xInterface-dxDE/2, xInterface+dxDE/2], ...
            [yInterface, yInterface], ...
            'k-', ...
            'LineWidth', 1.3);

    end
end

%% Sheet outline

hSheet = rectangle( ...
    'Position', [x0, y0, sheetLx, sheetLy], ...
    'EdgeColor', 'k', ...
    'LineWidth', 2);

%% Plot appearance

axis equal

xlim([x0, x0 + sheetLx]);
ylim([y0, y0 + sheetLy]);

xlabel('x [m]');
ylabel('y [m]');

title(sprintf( ...
    ['Grey: top 5 mm broken bonds | ', ...
     'Red: fully broken interfaces | ', ...
     'Black: floe boundaries'], ...
    nnz(actualFloeBoundary)), ...
    'FontWeight', 'bold');

set(gca, ...
    'FontSize', 14, ...
    'LineWidth', 1.2, ...
    'Box', 'on');

drawnow