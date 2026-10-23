clear;
clc;
close all;

%% =========================================================
% INPUT
% ==========================================================

resultFile = 'DE_floes_0600.mat';

% Minor ruler spacing in metres
minorTickStep_m = 0.10;


%% =========================================================
% LOAD FLOE RESULTS
% ==========================================================

if ~isfile(resultFile)
    error('Result file not found: %s', resultFile);
end

S = load(resultFile, 'Floes');
Floes = S.Floes;

requiredFields = { ...
    'ReferenceFile', ...
    'FloeID', ...
    'FloeTable', ...
    'DETable', ...
    'NumFloes'};

missingFields = setdiff(requiredFields, fieldnames(Floes));

if ~isempty(missingFields)
    error( ...
        'Floes structure is missing: %s', ...
        strjoin(missingFields, ', '));
end


%% =========================================================
% LOAD REFERENCE GEOMETRY
% ==========================================================

referenceFile = Floes.ReferenceFile;

% Try the result-file folder if the stored relative path is not found
if ~isfile(referenceFile)

    resultFolder = fileparts(resultFile);

    if isempty(resultFolder)
        resultFolder = pwd;
    end

    candidateFile = fullfile(resultFolder, referenceFile);

    if isfile(candidateFile)
        referenceFile = candidateFile;
    else
        error('Reference file not found: %s', Floes.ReferenceFile);
    end
end

R = load(referenceFile, 'Reference');
Reference = R.Reference;

numX = double(Reference.NumDE(1));
numY = double(Reference.NumDE(2));
numZ = double(Reference.NumDE(3));

dxDE = double(Reference.ElementSize(1));
dyDE = double(Reference.ElementSize(2));

xmin = double(Reference.Origin(1));
ymin = double(Reference.Origin(2));

if numZ ~= 1
    error('This top-view floe-size script requires NumDE(3) = 1.');
end

if dxDE <= 0 || dyDE <= 0
    error('The DE dimensions must be positive.');
end


%% =========================================================
% OUTPUT NAMES
% ==========================================================

[~, resultName] = fileparts(resultFile);

tagToken = regexp(resultName, '(\d+)$', 'tokens', 'once');

if isempty(tagToken)
    tagChar = resultName;
else
    tagChar = tagToken{1};
end

outDir = fullfile(pwd, sprintf('floe_outputs_%s', tagChar));

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

floeCsv = fullfile( ...
    outDir, ...
    sprintf('floe_statistics_%s.csv', tagChar));

afsdCsv = fullfile( ...
    outDir, ...
    sprintf('AFSD_%s.csv', tagChar));

statsFile = fullfile( ...
    outDir, ...
    sprintf('summary_%s.txt', tagChar));

afsdFigure = fullfile( ...
    outDir, ...
    sprintf('AFSD_%s.png', tagChar));

boundaryFigure = fullfile( ...
    outDir, ...
    sprintf('floe_boundaries_%s.png', tagChar));

scaleFigure = fullfile( ...
    outDir, ...
    sprintf('floe_boundaries_with_Leff_%s.png', tagChar));


%% =========================================================
% FLOE LABELS
% ==========================================================

floeID = double(Floes.FloeID(:));
numFloes = double(Floes.NumFloes);

nDE = numX*numY*numZ;

if numel(floeID) ~= nDE
    error( ...
        'Floes.FloeID contains %d values, but the grid has %d DEs.', ...
        numel(floeID), ...
        nDE);
end

if any(floeID < 1) || any(floeID > numFloes)
    error('Floes.FloeID contains invalid floe IDs.');
end

if max(floeID) ~= numFloes
    error( ...
        'Maximum FloeID is %d, but NumFloes is %d.', ...
        max(floeID), ...
        numFloes);
end

% DE IDs increase first in x and then in y.
% Transposing produces a numY-by-numX top-view grid.
labels = reshape(floeID, [numX, numY])';


%% =========================================================
% FLOE AREA
% ==========================================================

numDEPerFloe = accumarray( ...
    floeID, ...
    1, ...
    [numFloes, 1], ...
    @sum, ...
    0);

deArea = dxDE*dyDE;

floeAreas = numDEPerFloe*deArea;

totalArea = sum(floeAreas);

sheetArea = ...
    numX*dxDE * ...
    numY*dyDE;

areaDifference = totalArea - sheetArea;

if abs(areaDifference) > ...
        1e-10*max([totalArea, sheetArea, 1])

    warning( ...
        ['The summed floe area differs from the sheet area by ' ...
         '%.6g m^2.'], ...
        areaDifference);
end


%% =========================================================
% RECTANGLE-EQUIVALENT FLOE LENGTH
%
% W_f:
%   Complete y-span occupied by the floe.
%
% L_f:
%   Length of a rectangle having the same floe area and
%   the same y-span:
%
%       L_f = A_f / W_f
% ==========================================================

ySpanRows = zeros(numFloes, 1);
characteristicY_m = zeros(numFloes, 1);
floeLength_m = zeros(numFloes, 1);

for f = 1:numFloes

    floeMask = labels == f;

    % Rows containing at least one DE from this floe
    occupiedRows = find(any(floeMask, 2));

    if isempty(occupiedRows)
        error('Floe %d does not contain any DEs.', f);
    end

    % Full row span from the lowest occupied row to the highest
    numberOfRows = ...
        max(occupiedRows) - min(occupiedRows) + 1;

    ySpanRows(f) = numberOfRows;

    characteristicY_m(f) = ...
        numberOfRows*dyDE;

    floeLength_m(f) = ...
        floeAreas(f)/characteristicY_m(f);
end


%% =========================================================
% AREA WEIGHTS AND EFFECTIVE FLOE LENGTH
%
%       w_f = A_f / A_total
%
%       L_eff = sum(w_f L_f)
% ==========================================================

areaWeight = floeAreas/totalArea;

L_eff_m = sum( ...
    areaWeight .* floeLength_m);

weightSum = sum(areaWeight);


%% =========================================================
% AREA-BASED FLOE-SIZE DISTRIBUTION
%
% N = number of floes
%
% M = ceil(sqrt(N))
%
% delta_d = (d_max - d_min)/M
%
% AFSD = area fraction in bin / bin width
% ==========================================================

N = numFloes;
M = max(1, ceil(sqrt(N)));

dMin = min(floeLength_m);
dMax = max(floeLength_m);

sameLengthTolerance = ...
    eps(max([abs(dMin), abs(dMax), 1]));

if abs(dMax - dMin) <= sameLengthTolerance

    % All floes have the same equivalent length
    M = 1;

    % Give the single bin a finite physical width
    delta_d = dxDE;

    dEdges = [ ...
        dMin - delta_d/2, ...
        dMax + delta_d/2];

else

    delta_d = (dMax - dMin)/M;

    dEdges = ...
        dMin + (0:M)*delta_d;

    % Avoid floating-point exclusion of the maximum value
    dEdges(end) = dMax;
end

dCenters = ...
    0.5*(dEdges(1:end-1) + dEdges(2:end));

binWidths = diff(dEdges(:));

[binCount, ~, binID] = histcounts( ...
    floeLength_m, ...
    dEdges);

binCount = binCount(:);

validBin = binID > 0;

if nnz(validBin) ~= numFloes
    warning( ...
        '%d floes were not assigned to an AFSD bin.', ...
        numFloes - nnz(validBin));
end

binArea = accumarray( ...
    binID(validBin), ...
    floeAreas(validBin), ...
    [M, 1], ...
    @sum, ...
    0);

binAreaFraction = binArea/totalArea;

AFSD = binAreaFraction ./ binWidths;

AFSD_integral = sum(AFSD .* binWidths);


%% =========================================================
% SAVE FLOE STATISTICS
% ==========================================================

FloeStatistics = table( ...
    uint32((1:numFloes)'), ...
    uint32(numDEPerFloe), ...
    floeAreas, ...
    uint16(ySpanRows), ...
    characteristicY_m, ...
    floeLength_m, ...
    areaWeight, ...
    'VariableNames', { ...
        'FloeID', ...
        'NumDE', ...
        'Area_m2', ...
        'YSpanRows', ...
        'CharacteristicY_m', ...
        'EquivalentLength_m', ...
        'AreaWeight'});

writetable(FloeStatistics, floeCsv);


%% =========================================================
% SAVE AFSD DATA
% ==========================================================

AFSDTable = table( ...
    uint32((1:M)'), ...
    dEdges(1:end-1)', ...
    dEdges(2:end)', ...
    dCenters(:), ...
    binWidths, ...
    uint32(binCount), ...
    binArea, ...
    binAreaFraction, ...
    AFSD, ...
    'VariableNames', { ...
        'Bin', ...
        'LeftEdge_m', ...
        'RightEdge_m', ...
        'Center_m', ...
        'Width_m', ...
        'NumFloes', ...
        'Area_m2', ...
        'AreaFraction', ...
        'AFSD_m_inv'});

writetable(AFSDTable, afsdCsv);


%% =========================================================
% SAVE MATLAB RESULTS
% ==========================================================

Statistics.ResultFile = resultFile;
Statistics.ReferenceFile = referenceFile;

Statistics.FloeStatistics = FloeStatistics;
Statistics.AFSDTable = AFSDTable;

Statistics.TotalArea_m2 = totalArea;
Statistics.SheetArea_m2 = sheetArea;

Statistics.L_eff_m = L_eff_m;
Statistics.AFSDIntegral = AFSD_integral;

Statistics.NumFloes = numFloes;
Statistics.NumBins = M;
Statistics.DeltaD_m = delta_d;

statisticsMat = fullfile( ...
    outDir, ...
    sprintf('floe_statistics_%s.mat', tagChar));

save(statisticsMat, 'Statistics', '-v7.3');


%% =========================================================
% AFSD PLOT
% ==========================================================

fig1 = figure( ...
    'Color', 'w', ...
    'Position', [100, 100, 750, 500]);

bar(dCenters, AFSD, 1);

xlabel( ...
    'Rectangle-equivalent floe length [m]');

ylabel('P_A(d) [m^{-1}]');

title(sprintf( ...
    'AFSD: %d floes', ...
    numFloes), ...
    'FontWeight', 'bold');

grid on
box on

set(gca, ...
    'FontSize', 24, ...
    'LineWidth', 1.1);

exportgraphics( ...
    fig1, ...
    afsdFigure, ...
    'Resolution', 300);


%% =========================================================
% PLAIN FLOE MAP
% ==========================================================

fig2 = figure( ...
    'Color', 'w', ...
    'Position', [100, 100, 1300, 330]);

ax2 = axes('Parent', fig2);

drawFloeMap( ...
    ax2, ...
    labels, ...
    numX, ...
    numY, ...
    xmin, ...
    ymin, ...
    dxDE, ...
    dyDE, ...
    numFloes);

xlabel(ax2, 'x [m]');
ylabel(ax2, 'y [m]');

title(ax2, ...
    sprintf('Graph-derived floes: %d floes', numFloes), ...
    'FontWeight', 'bold');

set(ax2, ...
    'FontSize', 12, ...
    'LineWidth', 1.1);

exportgraphics( ...
    fig2, ...
    boundaryFigure, ...
    'Resolution', 300);


%% =========================================================
% FLOE MAP WITH RULER AND L_eff
% ==========================================================

fig3 = figure( ...
    'Color', 'w', ...
    'Position', [100, 100, 1500, 280]);

ax3 = axes( ...
    'Parent', fig3, ...
    'Position', [0.07, 0.18, 0.90, 0.72]);

drawFloeMapWithScale( ...
    ax3, ...
    labels, ...
    numX, ...
    numY, ...
    xmin, ...
    ymin, ...
    dxDE, ...
    dyDE, ...
    numFloes, ...
    L_eff_m, ...
    minorTickStep_m);

% Hide the normal y-axis scale
set(ax3, ...
    'YTick', [], ...
    'YTickLabel', [], ...
    'YColor', 'none', ...
    'FontSize', 12, ...
    'LineWidth', 1.1);

% Remove the outer axes box
box(ax3, 'off');% Hide the normal y-axis scale
set(ax3, ...
    'YTick', [], ...
    'YTickLabel', [], ...
    'YColor', 'none', ...
    'FontSize', 12, ...
    'LineWidth', 1.1);

% Remove the outer axes box
box(ax3, 'off');

exportgraphics( ...
    fig3, ...
    scaleFigure, ...
    'Resolution', 300);


%% =========================================================
% TEXT SUMMARY
% ==========================================================

fid = fopen(statsFile, 'w');

if fid < 0
    error('Could not create summary file: %s', statsFile);
end

fprintf(fid, 'FLOE STATISTICS SUMMARY\n');
fprintf(fid, '=======================\n\n');

fprintf(fid, 'Result file: %s\n', resultFile);
fprintf(fid, 'Reference file: %s\n\n', referenceFile);

fprintf(fid, 'Number of floes: %d\n', numFloes);
fprintf(fid, 'Number of DEs: %d\n', nDE);

fprintf(fid, 'DE size: %.8g x %.8g m\n', dxDE, dyDE);
fprintf(fid, 'Sheet size: %.8g x %.8g m\n', ...
    numX*dxDE, numY*dyDE);

fprintf(fid, 'Total floe area: %.12g m^2\n', totalArea);
fprintf(fid, 'Sheet area: %.12g m^2\n', sheetArea);
fprintf(fid, 'Area difference: %.12g m^2\n\n', areaDifference);

fprintf(fid, 'Area-weight sum: %.12g\n', weightSum);

fprintf(fid, 'L_eff: %.12g m\n', L_eff_m);
fprintf(fid, 'L_eff: %.12g cm\n\n', 100*L_eff_m);

fprintf(fid, 'AFSD bins: %d\n', M);
fprintf(fid, 'Nominal delta d: %.12g m\n', delta_d);
fprintf(fid, 'AFSD integral: %.12g\n', AFSD_integral);

fclose(fid);


%% =========================================================
% COMMAND-WINDOW SUMMARY
% ==========================================================

fprintf('\nFLOE-SIZE RESULTS\n');
fprintf('Number of floes: %d\n', numFloes);
fprintf('Total area: %.6g m^2\n', totalArea);
fprintf('Sheet area: %.6g m^2\n', sheetArea);

fprintf('\nRectangle-equivalent lengths\n');
fprintf('Minimum length: %.6g m\n', dMin);
fprintf('Maximum length: %.6g m\n', dMax);
fprintf('L_eff: %.6g m\n', L_eff_m);
fprintf('L_eff: %.6g cm\n', 100*L_eff_m);

fprintf('\nAFSD\n');
fprintf('Number of bins: M = %d\n', M);
fprintf('Delta d: %.6g m\n', delta_d);
fprintf('AFSD integral: %.8f\n', AFSD_integral);


%% =========================================================
% LOCAL FUNCTIONS
% ==========================================================

function drawFloeMap( ...
    ax, labels, Nx, Ny, xmin, ymin, ...
    dxDE, dyDE, numFloes)

    hold(ax, 'on');
    axis(ax, 'equal');
    box(ax, 'on');

    colormap(ax, lines(max(numFloes, 1)));

    if numFloes == 1
        clim(ax, [0.5, 1.5]);
    else
        clim(ax, [0.5, numFloes + 0.5]);
    end

    %% Color every DE according to its floe ID
    for j = 1:Ny

        for i = 1:Nx

            floeNumber = labels(j, i);

            x1 = xmin + (i - 1)*dxDE;
            x2 = xmin + i*dxDE;

            y1 = ymin + (j - 1)*dyDE;
            y2 = ymin + j*dyDE;

            patch( ...
                ax, ...
                [x1, x2, x2, x1], ...
                [y1, y1, y2, y2], ...
                floeNumber, ...
                'EdgeColor', 'none');
        end
    end

    %% Vertical boundaries between different floes
    for j = 1:Ny

        for i = 1:Nx-1

            if labels(j, i) ~= labels(j, i+1)

                xBoundary = xmin + i*dxDE;

                y1 = ymin + (j - 1)*dyDE;
                y2 = ymin + j*dyDE;

                plot( ...
                    ax, ...
                    [xBoundary, xBoundary], ...
                    [y1, y2], ...
                    'k-', ...
                    'LineWidth', 1.8);
            end
        end
    end

    %% Horizontal boundaries between different floes
    for j = 1:Ny-1

        for i = 1:Nx

            if labels(j, i) ~= labels(j+1, i)

                yBoundary = ymin + j*dyDE;

                x1 = xmin + (i - 1)*dxDE;
                x2 = xmin + i*dxDE;

                plot( ...
                    ax, ...
                    [x1, x2], ...
                    [yBoundary, yBoundary], ...
                    'k-', ...
                    'LineWidth', 1.8);
            end
        end
    end

    %% Outer sheet boundary
    sheetWidth = Nx*dxDE;
    sheetHeight = Ny*dyDE;

    rectangle( ...
        ax, ...
        'Position', [ ...
            xmin, ...
            ymin, ...
            sheetWidth, ...
            sheetHeight], ...
        'EdgeColor', 'k', ...
        'LineWidth', 2.2);

    xlim(ax, [xmin, xmin + sheetWidth]);
    ylim(ax, [ymin, ymin + sheetHeight]);

    set(ax, 'YDir', 'normal');
end


function drawFloeMapWithScale( ...
    ax, labels, Nx, Ny, xmin, ymin, ...
    dxDE, dyDE, numFloes, L_eff_m, minorStep_m)

    %% Draw colored floes and boundaries
    drawFloeMap( ...
        ax, ...
        labels, ...
        Nx, ...
        Ny, ...
        xmin, ...
        ymin, ...
        dxDE, ...
        dyDE, ...
        numFloes);

    % drawFloeMap uses equal x-y scaling.
    % Override it for this wide ruler figure.
    axis(ax, 'normal');
    daspect(ax, 'auto');

    plateLength = Nx*dxDE;
    plateHeight = Ny*dyDE;

    xLeft = xmin;
    xRight = xmin + plateLength;

    yBottom = ymin;
    yTop = ymin + plateHeight;


    %% Major ruler ticks
    majorStep = niceStep(plateLength/9);

    majorDistances = 0:majorStep:plateLength;

    if isempty(majorDistances) || ...
            abs(majorDistances(end) - plateLength) > 1e-12

        majorDistances = [majorDistances, plateLength];
    end

    majorTicks = xLeft + majorDistances;


    %% Minor ruler ticks
    minorDistances = 0:minorStep_m:plateLength;

    if isempty(minorDistances) || ...
            abs(minorDistances(end) - plateLength) > 1e-12

        minorDistances = [minorDistances, plateLength];
    end

    minorTicks = xLeft + minorDistances;


    %% Ruler dimensions
    minorLength = 0.08*plateHeight;
    majorLength = 0.14*plateHeight;
    labelGap = 0.035*plateHeight;

    rulerTextY = ...
        yBottom - majorLength - labelGap;


    %% Minor ticks
    for xTick = minorTicks

        plot( ...
            ax, ...
            [xTick, xTick], ...
            [yBottom, yBottom - minorLength], ...
            'k-', ...
            'LineWidth', 0.6);
    end


    %% Major ticks and labels
    for k = 1:numel(majorTicks)

        xTick = majorTicks(k);
        distanceFromLeft = majorDistances(k);

        plot( ...
            ax, ...
            [xTick, xTick], ...
            [yBottom, yBottom - majorLength], ...
            'k-', ...
            'LineWidth', 1.0);

        text( ...
            ax, ...
            xTick, ...
            rulerTextY, ...
            sprintf('%.3g', distanceFromLeft), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'top',...
            'FontSize',24);
    end


    %% L_eff marker above the sheet
    leffY = yTop + 0.40*plateHeight;

    leffTickLength = 0.08*plateHeight;
    leffTextY = leffY + 0.035*plateHeight;

    xLeff = min(xLeft + L_eff_m, xRight);

    plot( ...
        ax, ...
        [xLeft, xLeff], ...
        [leffY, leffY], ...
        'k-', ...
        'LineWidth', 3.0);

    plot( ...
        ax, ...
        [xLeft, xLeft], ...
        [ ...
            leffY - leffTickLength/2, ...
            leffY + leffTickLength/2], ...
        'k-', ...
        'LineWidth', 1.3);

    plot( ...
        ax, ...
        [xLeff, xLeff], ...
        [ ...
            leffY - leffTickLength/2, ...
            leffY + leffTickLength/2], ...
        'k-', ...
        'LineWidth', 1.3);

    text( ...
        ax, ...
        1.2*(xLeft + xLeff), ...
        leffTextY, ...
        sprintf( ...
            'L_{eff} = %.4g m (%.4g cm)', ...
            L_eff_m, ...
            100*L_eff_m), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontWeight', 'bold',...
        'FontSize',18);


    %% Compact vertical limits
    extraBottom = 0.7*plateHeight;
    extraTop = 0.50*plateHeight;

    xlim(ax, [xLeft, xRight]);

    ylim(ax, [ ...
        yBottom - extraBottom, ...
        yTop + extraTop]);

    % Custom ruler replaces the normal x-axis ticks
    set( ...
        ax, ...
        'XTick', [], ...
        'XTickLabel', [], ...
        'YDir', 'normal',...
        'FontSize',12);
    axis(ax, 'normal');
    % daspect(ax, 'auto');
    box(ax, 'off');
end


function step = niceStep(rawStep)

    if rawStep <= 0 || ~isfinite(rawStep)
        step = 1;
        return
    end

    scale = 10^floor(log10(rawStep));

    choices = [1, 2, 5, 10]*scale;

    index = find( ...
        choices >= rawStep, ...
        1, ...
        'first');

    if isempty(index)
        step = 10*scale;
    else
        step = choices(index);
    end
end