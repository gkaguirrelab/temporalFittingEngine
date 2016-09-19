% retSearch.m

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

%convert degrees to pixels
stimSize = 39.2257;
subjectDist = 106.5;
screenHgt = 1080;
DVA = rad2deg(2*atan(stimSize/(2*subjectDist)));
pxlPerDeg = round(((screenHgt)/(DVA))/2);

X = 0:pxlPerDeg:size(downStim,1);
Y = 0:pxlPerDeg:size(downStim,2);
X = X(2:end-1);
Y = Y(2:end-1);
sigmaList = 13:13:26*5;



%% generat a predictions
count = 1;
tic 
for x = 1:length(X)
    for y = 1:length(Y)
        for s = 1:length(sigmaList)

            TCcell.TCmat(count,:) = makePredTC(downStim,X(x),Y(y),sigmaList(s));
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









