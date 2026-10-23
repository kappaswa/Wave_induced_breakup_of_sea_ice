
%% ================================================================

clear; clc; close all;

%% ---------------- User settings ----------------

folder = './';
filePattern = '2d_*.csv';

dp = 0.005;                  % particle spacing [m]
binWidth = 5*dp;             % x-bin width [m]
minPtsBin = 2;               % minimum particles in bin

% Surface extraction
surfaceMethod = 'percentile';    % 'max' or 'percentile'
surfacePercentile = 95;          % safer than max for noisy SPH surface

% Overwash removal
tolBottom = 0.5*dp;          % only use water below float bottom - tolBottom

% Crack handling
includeCrackWaterAsSurface = false;
crackGapBuffer = 4*dp;       % used to decide ice-covered x-region

% Gap filling / smoothing
maxGapToFill = 6;            % fill only small NaN gaps, in bins
smoothSurface = true;
smoothBins = 3;

% Movie settings
makeVideo = true;
videoName = 'eta_no_overwash_movie.mp4';
videoFPS = 12;

% CSV output
outCSV = 'eta_no_overwash_vs_x_time.csv';

% Plot limits
xlimManual = [];
ylimManual = [];             % example: [-0.05 0.07]
ylimPad = 0.005;             % padding around eta curve if ylimManual = []

%% ---------------- Find and sort files ----------------

files = dir(fullfile(folder, filePattern));

if isempty(files)
    error('No files found matching: %s', fullfile(folder, filePattern));
end

fileNames = {files.name};
fileNums = nan(size(fileNames));

for i = 1:numel(fileNames)
    tokens = regexp(fileNames{i}, '(\d+)(?=\.csv$)', 'tokens');
    if ~isempty(tokens)
        fileNums(i) = str2double(tokens{end}{1});
    else
        fileNums(i) = i;
    end
end

[~, order] = sort(fileNums);
files = files(order);
nFiles = numel(files);

fprintf('Found %d files.\n', nFiles);

%% ---------------- First pass: common x-grid ----------------

xminAll = inf;
xmaxAll = -inf;

fprintf('First pass: finding x limits...\n');

for k = 1:nFiles
    fname = fullfile(files(k).folder, files(k).name);

    A = readmatrix(fname, ...
        'FileType','text', ...
        'Delimiter',';', ...
        'NumHeaderLines',4);

    A = A(:,1:10);

    x = A(:,1);
    type = A(:,10);

    isFluid = type == 3;
    isFloat = type == 2;

    xUse = x(isFluid | isFloat);

    xminAll = min(xminAll, min(xUse));
    xmaxAll = max(xmaxAll, max(xUse));
end

edges = xminAll:binWidth:(xmaxAll + binWidth);
xGrid = edges(1:end-1) + binWidth/2;
nx = numel(xGrid);

if isempty(xlimManual)
    xlimUse = [min(xGrid), max(xGrid)];
else
    xlimUse = xlimManual;
end

%% ---------------- Storage ----------------

time = nan(nFiles,1);
zSurfaceAll = nan(nFiles,nx);
etaAll = nan(nFiles,nx);

%% ---------------- Main extraction loop ----------------

fprintf('Extracting overwash-free surface...\n');

for k = 1:nFiles

    fname = fullfile(files(k).folder, files(k).name);

    time(k) = readTimeFromDualSPHysicsCSV(fname);

    A = readmatrix(fname, ...
        'FileType','text', ...
        'Delimiter',';', ...
        'NumHeaderLines',4);

    A = A(:,1:10);

    x = A(:,1);
    z = A(:,3);
    type = A(:,10);

    isFluid = type == 3;
    isFloat = type == 2;

    xf = x(isFluid);
    zf = z(isFluid);

    xb = x(isFloat);
    zb = z(isFloat);

    zSurface = getSurfaceNoOverwashSimple( ...
        xf, zf, xb, zb, edges, ...
        minPtsBin, tolBottom, ...
        surfaceMethod, surfacePercentile, ...
        includeCrackWaterAsSurface, crackGapBuffer, ...
        maxGapToFill);

    if smoothSurface
        zSurface = smoothContinuousSegments(zSurface, smoothBins);
    end

    zSurfaceAll(k,:) = zSurface;

    fprintf('Processed %d / %d: %s, t = %.4f s\n', ...
        k, nFiles, files(k).name, time(k));
end

%% ---------------- Elevation relative to initial ----------------

zInitial = zSurfaceAll(1,:);

for k = 1:nFiles
    etaAll(k,:) = zSurfaceAll(k,:) - zInitial;
end

%% ---------------- Write CSV ----------------
% Format:
%   row 1: NaN, x1, x2, x3, ...
%   row 2: t1, eta(x1,t1), eta(x2,t1), ...
%   row 3: t2, eta(x1,t2), eta(x2,t2), ...

outMatrix = [nan, xGrid; time, etaAll];

writematrix(outMatrix, outCSV);

fprintf('Saved CSV: %s\n', outCSV);

%% ---------------- Make simple elevation movie ----------------

if makeVideo

    etaFinite = etaAll(isfinite(etaAll));

    if isempty(ylimManual)
        if isempty(etaFinite)
            yLimits = [-0.05 0.05];
        else
            ymin = min(etaFinite);
            ymax = max(etaFinite);

            pad = max(ylimPad, 0.10*(ymax - ymin + eps));
            yLimits = [ymin - pad, ymax + pad];
        end
    else
        yLimits = ylimManual;
    end

    vw = VideoWriter(videoName, 'Motion JPEG AVI');
    vw.FrameRate = videoFPS;
    open(vw);

    fig = figure('Color','w', 'Position', [100 100 1100 650]);

    for k = 1:nFiles

        clf(fig);
        hold on; box on; grid on;

        plot(xGrid, etaAll(k,:), 'o-', ...
            'LineWidth', 1.5, ...
            'MarkerSize', 5);

        yline(0, 'k--');

        xlabel('x [m]');
        ylabel('\eta = elevation - initial elevation [m]');
        title(sprintf('Elevation against x, Part = %d', k));

        if isnan(time(k))
            labelText = sprintf('Part = %d', k);
        else
            labelText = sprintf('t = %.3f s', time(k));
        end

        text(0.03, 0.90, labelText, ...
            'Units','normalized', ...
            'FontSize',16, ...
            'FontWeight','bold', ...
            'BackgroundColor','w', ...
            'EdgeColor','k');

        xlim(xlimUse);
        ylim(yLimits);

        drawnow;

        frame = getframe(fig);
        writeVideo(vw, frame);
    end

    close(vw);

    fprintf('Saved movie: %s\n', videoName);
end

fprintf('Done.\n');

%% ========================================================================
%  Function: get simple water surface excluding overwash
%% ========================================================================

function zSurface = getSurfaceNoOverwashSimple( ...
    xf, zf, xb, zb, edges, ...
    minPtsBin, tolBottom, ...
    surfaceMethod, surfacePercentile, ...
    includeCrackWaterAsSurface, crackGapBuffer, ...
    maxGapToFill)

    nb = numel(edges) - 1;
    xCenters = edges(1:end-1) + diff(edges)/2;

    zSurface = nan(1,nb);
    zFloatBottom = nan(1,nb);
    hasFloat = false(1,nb);

    %% Overall ice-covered region
    if ~isempty(xb)
        iceXmin = min(xb) - crackGapBuffer;
        iceXmax = max(xb) + crackGapBuffer;
        insideIceRegion = xCenters >= iceXmin & xCenters <= iceXmax;
    else
        insideIceRegion = false(size(xCenters));
    end

    %% Bin-by-bin surface extraction
    for i = 1:nb

        Ifloat = xb >= edges(i) & xb < edges(i+1);
        Ifluid = xf >= edges(i) & xf < edges(i+1);

        if nnz(Ifloat) >= minPtsBin
            zFloatBottom(i) = min(zb(Ifloat));
            hasFloat(i) = true;
        end

        if nnz(Ifluid) < minPtsBin
            continue;
        end

        zLocalFluid = zf(Ifluid);

        if hasFloat(i)
            % Main overwash removal:
            % do not use fluid on top of the float.
            keep = zLocalFluid <= zFloatBottom(i) - tolBottom;

            if nnz(keep) >= minPtsBin
                zSurface(i) = surfaceValue(zLocalFluid(keep), surfaceMethod, surfacePercentile);
            else
                zSurface(i) = nan;
            end

        else
            % No float in this x-bin.
            % Inside ice region this is likely a crack/gap, so do not let
            % the waveform climb vertically unless user asks for that.
            if insideIceRegion(i) && ~includeCrackWaterAsSurface
                zSurface(i) = nan;
            else
                zSurface(i) = surfaceValue(zLocalFluid, surfaceMethod, surfacePercentile);
            end
        end
    end

    %% Fill only small gaps
    zSurface = fillSmallNaNGaps(zSurface, xCenters, maxGapToFill);
end

%% ========================================================================
%  Function: choose max or percentile surface value
%% ========================================================================

function val = surfaceValue(zLocal, surfaceMethod, surfacePercentile)

    switch lower(surfaceMethod)
        case 'max'
            val = max(zLocal);

        case 'percentile'
            val = prctile(zLocal, surfacePercentile);

        otherwise
            error('surfaceMethod must be max or percentile.');
    end
end

%% ========================================================================
%  Function: fill only small NaN gaps
%% ========================================================================

function y = fillSmallNaNGaps(y, x, maxGapToFill)

    nanRun = isnan(y);

    d = diff([false, nanRun, false]);
    runStart = find(d == 1);
    runEnd = find(d == -1) - 1;
    runLength = runEnd - runStart + 1;

    for r = 1:numel(runStart)

        if runLength(r) <= maxGapToFill

            idx = runStart(r):runEnd(r);

            leftIdx = runStart(r) - 1;
            rightIdx = runEnd(r) + 1;

            if leftIdx >= 1 && rightIdx <= numel(y) && ...
                    isfinite(y(leftIdx)) && isfinite(y(rightIdx))

                y(idx) = interp1( ...
                    x([leftIdx rightIdx]), ...
                    y([leftIdx rightIdx]), ...
                    x(idx), ...
                    'linear');
            end
        end
    end
end

%% ========================================================================
%  Function: smooth only continuous non-NaN segments
%% ========================================================================

function ySmooth = smoothContinuousSegments(y, smoothBins)

    ySmooth = y;

    if smoothBins <= 1
        return;
    end

    valid = isfinite(y);

    d = diff([false, valid, false]);
    segStart = find(d == 1);
    segEnd = find(d == -1) - 1;

    for s = 1:numel(segStart)

        idx = segStart(s):segEnd(s);

        if numel(idx) >= smoothBins
            ySmooth(idx) = movmean(y(idx), smoothBins, 'omitnan');
        end
    end
end

%% ========================================================================
%  Function: read timestep from DualSPHysics CSV header
%% ========================================================================

function t = readTimeFromDualSPHysicsCSV(fname)

    t = nan;

    fid = fopen(fname, 'r');

    if fid < 0
        return;
    end

    try
        fgetl(fid);          % header line
        line2 = fgetl(fid);  % values line
        fclose(fid);

        parts = split(string(line2), ';');

        if numel(parts) >= 1
            t = str2double(parts(1));
        end

    catch
        fclose(fid);
        t = nan;
    end
end