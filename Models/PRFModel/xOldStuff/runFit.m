% tests the predictions for the retinotopy model

%% set defaults
dataDir = '/Users/Shared/Matlab/gkaguirrelab/mriTemporalFitting/PRFToolbox/';

%% load the stimulus file
tmp = load(fullfile(dataDir,'pRFimages.mat'));

%% Binarize the stimulus
stim = tmp.imagesFull;
stim(stim ~=128) = 1;
stim(stim == 128 ) = 0;

%% Downsample the frames 
framesPerPos = 8;
start = 1:framesPerPos:size(stim,3);
stop = start(2)-1:framesPerPos:size(stim,3);
for i = 1:length(start)
    downStim(:,:,i) = mean(stim(:,:,start(i):stop(i)),3);
end

%% Downsample stim size 

%convert degrees to pixels
stimSize = 39.2257;
subjectDist = 106.5;
screenHgt = 1080;
DVA = rad2deg(2*atan(stimSize/(2*subjectDist)));
pxlPerDeg = ((screenHgt)/(DVA))/2;
resample2deg = round(pxlPerDeg/0.5);

smallStim = imresize(downStim,[resample2deg NaN]);
smallStim(smallStim >=0.4) = 1;
smallStim(smallStim<0.4) = 0;


%% Make simulated timecourse 
%TC = makePredTC(smallStim,25,60,1);
load V1tc.mat;
tic

for voxel = 1:size(V1tc,1)

    TC = V1tc(voxel,:);
    %% Grid search for starting params
 
    X = 0:floor(size(smallStim,1)/3):size(smallStim,1);
    Y = 0:floor(size(smallStim,2)/3):size(smallStim,2);
    X = X(2:end-1);
    Y = Y(2:end-1);
    for x = 1:length(X)
        for y = 1:length(Y);
            paramsVec=[X(x),Y(y),5];
            corrMat(x,y) = retCorr(smallStim,TC,paramsVec);
        end
    end

    [a,b] = find(corrMat == min(corrMat(:)));
    initX = X(a);
    initY = Y(b);


    %% test the fit
    R = @(paramsVec)retCorr(smallStim,TC,paramsVec);
    preds(voxel,:) =fmincon(R,[initX, initY, 5],[],[],[],[],[0 0 .4],[size(smallStim,1)*1.5,size(smallStim,2)*1.5,30]);

end
toc
finalPreds = preds/2;
[th ecc] = cart2pol(finalPreds(:,1),finalPreds(:,2));
th = rad2deg(th);
results = [ecc,th,finalPreds(:,3)];



save('results', 'results', '-v7.3')