% retSearch.m
%
%   Script to find pRFs
%
%   Written by Michael Barnett and Andrew S Bock Sep 2016

%% set defaults
[~, tmpName]            = system('whoami');
userName                = strtrim(tmpName); % Get user name
dataDir                 = ['/Users/' userName '/Dropbox-Aguirre-Brainard-Lab/retData'];
stimSize                = 39.2257;
fieldSize               = stimSize/2; % radius of stimuluated visual field in degrees visual angle
subjectDist             = 106.5;
screenHgt               = 1080;
framesPerTR             = 8;
gridScale               = 1;
maxXY                   = fieldSize*gridScale;
gridPoints              = 50; % based on mrVista (default is 50)
sampleRate              = maxXY./gridPoints; % sample rate in visual angle
sigList                 = 1;
TR                      = 0.8; % seconds
%% load the stimulus file
tmp                     = load(fullfile(dataDir,'pRFimages.mat'));

%% Binarize the stimulus
stim                    = tmp.imagesFull;
stim(stim ~=128)        = 1;
stim(stim == 128 )      = 0;

%% Average the frames within each TR
start                   = 1:framesPerTR:size(stim,3);
stop                    = start(2)-1:framesPerTR:size(stim,3);
meanImages              = nan(size(stim,1),size(stim,2),size(stim,3)/framesPerTR);
for i = 1:length(start)
    meanImages(:,:,i) = mean(stim(:,:,start(i):stop(i)),3);
end
%% Create X, Y, and sigma
tmpgrid                 = -maxXY:sampleRate:maxXY;
[x,y]                   = meshgrid(tmpgrid,tmpgrid);
tmpx0                   = x(:);
tmpy0                   = y(:);
X                       = x(:);
Y                       = y(:);
x0                      = repmat(tmpx0,size(sigList,1),1);
y0                      = repmat(tmpy0,size(sigList,1),1);
sigs                    = repmat(sigList,size(tmpx0,1),1);
%% resample images to sampling grid
nImages = size(meanImages, 3);
resampled = zeros((1+gridScale*(gridPoints/gridScale))^2,nImages);
for ii = 1:nImages
    tmp_im = imresize(meanImages(:,:,ii), 1+gridScale*[gridPoints/gridScale gridPoints/gridScale]);
    resampled(:, ii) = tmp_im(:);
end
%% Add black around stimulus region, to model the actual visual field (not just the bars)

%%% need to do this %%%

images = resampled; % do this for now
%% Break up into smaller matrices
nn = numel(x0); % grid points
[predPerTask,predTasks] = calc_tasks(nn,ceil(nn/1000));
predidx = [];
for i = 1:predTasks
    if isempty(predidx);
        predidx = [1,predPerTask(i)];
    else
        predidx = [predidx;[predidx(end,2)+1,predidx(end,2)+predPerTask(i)]];
    end
    predvals{i} = predidx(i,1):predidx(i,2);
end
%% Make/load HRF
HRF = doubleGammaHrf(TR);
%% Make predicted timecoures from stimulus images
for n=1:length(predvals)
    nSigs               = sigs(predvals{n},:);
    nx0                 = x0(predvals{n});
    ny0                 = y0(predvals{n});
    % Allow x, y, sigma to be a matrix so that the final output will be
    % size(X,1) by size(x0,2). This way we can make many RFs at the same time.
    if numel(nSigs)~=1,
        sz1 = size(X,1);
        sz2 = size(nSigs,1);
        X   = repmat(X,1,sz2);
        Y   = repmat(Y(:),1,sz2);
        nx0 = repmat(nx0',sz1,1);
        ny0 = repmat(ny0',sz1,1);
        nSigs = repmat(nSigs,1,1,sz1);
        nSigs = permute(nSigs,[3 1 2]);
    end
    % Translate grid so that center is at RF center
    nX = X - nx0;   % positive x0 moves center right
    nY = Y - ny0;   % positive y0 moves center up
    % make gaussian on current grid
    rf = exp (-(nY.^2 + nX.^2) ./ (2*nSigs(:,:,1).^2));
    % Convolve images with HRF
    allstimimages = filter(HRF,1, images');
    % Convolve images (with HRF) with receptive field
    pred = allstimimages*rf;
    % Set timecourses with very little variation (var<0.1) to flat
    % above images are set to be %BOLD/degree2
    pred = set_to_flat(pred);
    % store
    predTCs(:,(predvals{n} + (h-1)*nn)) = pred;
end
%%





%% Deprecated



%% convert degrees to pixels
DVA = rad2deg(2*atan(stimSize/(2*subjectDist)));
pxlPerDeg = round(((screenHgt)/(DVA))/2);
X = 0:pxlPerDeg:size(images,1);
Y = 0:pxlPerDeg:size(images,2);
X = X(2:end-1);
Y = Y(2:end-1);
%% generat a predictions
count = 1;
tic
for x = 1:length(X)
    for y = 1:length(Y)
        for s = 1:length(sigmaList)
            
            TCcell.TCmat(count,:) = makePredTC(images,X(x),Y(y),sigmaList(s));
            TCcell.params(count,:) = [X(x),Y(y),sigmaList(s)];
            count = count+1;
        end
    end
end
toc
save('TCcell', 'TCcell', '-v7.3')
% V1tc

for v = 1:size(V1tc,1)
    r = corr(V1tc(v,:)',TCcell.TCmat');
    paramPreds(v,:) = TCcell.params(find( r == max(r)),:);
end









