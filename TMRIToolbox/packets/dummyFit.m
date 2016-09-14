function [B,R2]         = dummyFit(packet,eventNum)

% Example function to obtain the beta and r-squared values for a given event

%% Convolve the stimulus event with an HRF
predResp                = filter(packet.kernel.values,1,packet.stimulus.values(eventNum,:));
% Downsample to the resolution of the data
nSamps                  = size(packet.stimulus.values,2) / size(packet.response.values,2);
downPredResp            = downsample(predResp,nSamps);
%% Run linear regression
X                       = [ones(size(downPredResp));downPredResp]';
Y                       = packet.response.values';
[B,~,~,~,STATS]         = regress(Y,X);
B                       = B(2);
R2                      = STATS(1);