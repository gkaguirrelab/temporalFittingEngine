function [fVal, responsePredicted] = crossValBDCM(obj,thePackets,leaveOutInd)

% function [fVal, responsePredicted] = crossValBDCM(obj,thePackets,leaveOutInd)
%
% implements BDCM cross-validation. Assumes same HRF used for all packets,
% same stimulus timebase, same response timebase, etc.
%
% inputs
%
% thePackets : list of packets
% leaveOutInd: indices of packets to leave out for fitting, and to use for
%              cross-validation
%
% outputs
%
% fVal             : error amount
% responsePredicted: fitted response

p = inputParser;
p.addRequired('thePackets',@iscell);
p.addRequired('leaveOutInd',@isnumeric);
p.parse(thePackets,leaveOutInd);

% get packet indices
packetInds = 1:length(thePackets);
% grab packets to fit
packetsToFitInd = packetInds(packetInds~=leaveOutInd);
packetsToFit = thePackets(packetsToFitInd);

% concatenate packets
packetsConcToFit = obj.concatenatePackets(packetsToFit);
% put HRF back in
packetsConcToFit.HRF = thePackets{1}.HRF;
% put in stimulus metadata
metaDataConcFit = obj.concMetaDataBDCM(packetsToFit);
packetsConcToFit.stimulus.metaData.theFrequencyIndices = metaDataConcFit.theFrequencyIndices;
% locking matrix
paramLockMatrix = createParamLockMatrixVanilla(unique(packetsConcToFit.stimulus.metaData.theFrequencyIndices) ...
                                               ,packetsConcToFit.stimulus.metaData.theFrequencyIndices,2);

defaultParamsInfo.nEvents = size(packetsConcToFit.stimulus.values,1);

%% Set parameters

params0 = obj.defaultParams('DefaultParamsInfo',defaultParamsInfo);
fprintf('Default model parameters:\n');
obj.print(params0);

%% Test paramsToVec and vecToParams
params1 = params0;
params1.paramsMainMatrix = rand(size(params1.paramMainMatrix));
x1 = obj.paramsToVec(params1);
params2 = obj.vecToParams(x1);
if (any(params1.paramMainMatrix ~= params2.paramMainMatrix))
    error('vecToParams and paramsToVec do not invert');
end

% do the fit
[paramsFit,fVal,fitResponse] = obj.fitResponse(packetsConcToFit,'DefaultParamsInfo',defaultParamsInfo, ...
                          'paramLockMatrix',paramLockMatrix);
fprintf('Model parameter from fits:\n');

% get mean parameter values for plotting
[~, meanParamValues,~] = ...
obj.plotParams(paramsFit,packetsConcToFit.stimulus.metaData.theFrequencyIndices,0);

% plot amplitudes
theFreq = thePackets{1}.stimulus.metaData.params.theFrequenciesHz ...
           (2:length(thePackets{1}.stimulus.metaData.params.theFrequenciesHz));
theFreqForPlot = [1 theFreq];
figure;
plot(theFreqForPlot,meanParamValues(1,:),'-kd','LineWidth',2,'MarkerSize',15); hold on
set(gca,'Xscale','log'); xlabel('Temporal Frequency (Hz)'); ylabel('Amplitude');
set(gca,'FontSize',15); set(gca,'Xtick',theFreq); axis square;
%% CROSS-VALIDATION
% sort parameters for cross-validation
paramsSorted = paramSort4crossVal(obj,meanParamValues, ...
    unique(packetsConcToFit.stimulus.metaData.theFrequencyIndices), ...
           packetsConcToFit.stimulus.metaData.theFrequencyIndices');
       
% grab packets for cross-validation 
packetsToCrossVal = thePackets(leaveOutInd);
% concatenate
packetsConcToCrossVal = obj.concatenatePackets(packetsToCrossVal);
% compute cross-validated fit
[fVal,responsePredicted] = fitError(obj,paramsSorted, ...
    packetsConcToCrossVal.stimulus,packetsConcToCrossVal.response,thePackets{1}.HRF);

obj.plot(packetsConcToCrossVal.response.timebase,packetsConcToCrossVal.response.values); hold on;
obj.plot(packetsConcToCrossVal.response.timebase,responsePredicted,'Color',[0 1 0],'NewWindow',false);

end