% % retSearch.m
% 
% %% set defaults
% dataDir = '/Users/Shared/Matlab/gkaguirrelab/mriTemporalFitting/PRFToolbox/';
% 
% %% load the stimulus file
% tmp = load(fullfile(dataDir,'pRFimages.mat'));

% %% Binarize the stimulus
% stim = tmp.imagesFull;
% stim(stim ~=128) = 1;
% stim(stim == 128 ) = 0;

%load binary stim images

load('/Users/micalan/Dropbox (Aguirre-Brainard Lab)/retData/binImageFull.mat');

%% Downsample the frames 
framesPerPos = 8;
start = 1:framesPerPos:size(binImageFull,3);
stop = start(2)-1:framesPerPos:size(binImageFull,3);
for i = 1:length(start)
    downStim(:,:,i) = mean(binImageFull(:,:,start(i):stop(i)),3);
end
padStim = padarray(downStim,[round(0.5*size(downStim,1)),round(0.5*size(downStim,1))]);
clear downStim binImageFull

%convert degrees to pixels
stimSize = 39.2257;
subjectDist = 106.5;
screenHgt = 1080;
DVA = rad2deg(2*atan(stimSize/(2*subjectDist)));
pxlPerDeg = round(((screenHgt)/(DVA))/2);

X = 0:2*pxlPerDeg:size(padStim,1);
Y = 0:2*pxlPerDeg:size(padStim,2);
X = X(2:end-1);
Y = Y(2:end-1);
sigmaList = [26];



%% generate a predictions
count = 1;
tic 
for x = 1:length(X)
    display(x)
    for y = 1:length(Y)
        display(y)
        for s = 1:length(sigmaList)

            TCcell.TCmat(count,:) = makePredTC(padStim,X(x),Y(y),sigmaList(s));
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









