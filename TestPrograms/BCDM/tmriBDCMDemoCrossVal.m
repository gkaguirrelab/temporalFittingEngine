%% LOAD PACKETS

load('/Users/benchin/Desktop/MELA_data/packets_asb1_041416.mat');
%% TRY CROSS-VALIDATION FOR JUST THE L-M
indicesOfInterest = [2 5 8 11];

thePackets = {};

for i = 1:length(indicesOfInterest)
   aPacket = packets{indicesOfInterest(i)};
   aPacket.response.values = aPacket.response.values{indicesOfInterest(i)};
   thePackets{length(thePackets)+1} = aPacket; 
end

clear packets;

%% Construct the model object
tmri = tmriBDCM;

%% DOWNSAMPLE STIMULI TO THE DESIRED RESOLUTION

% common factor to downsample by
timeFactor = 200;

startTime = thePackets{1}.stimulus.timebase(1);
endTime = thePackets{1}.stimulus.timebase(end);
nTimeSamples = round(endTime/timeFactor);
newStimulusTimebase = linspace(startTime,endTime,nTimeSamples+1);

%% DOWNSAMPLE HRF TO THE DESIRED RESOLUTION

startTime = thePackets{1}.HRF.timebase(1);
endTime = thePackets{1}.HRF.timebase(end);
nTimeSamples = round(endTime/timeFactor);
newHRFTimebase = linspace(startTime,endTime,nTimeSamples+1);

%%
% May want a utility function to convert lists of packets into lists of
% their substructures, and back.
for ii = 1:length(thePackets)
    theStimulusList{ii} = thePackets{ii}.stimulus;
    theHRFList{ii} = thePackets{ii}.HRF;
end

resampledStimulusList = tmri.resamplePacketStruct(theStimulusList,newStimulusTimebase);
resampledHRFList = tmri.resamplePacketStruct(theHRFList,newHRFTimebase);

for ii = 1:length(thePackets)
    thePackets{ii}.stimulus = resampledStimulusList{ii};
    theHRFList{ii} = resampledHRFList{ii};
    thePackets{ii}.HRF = theHRFList{ii};
end

%%
% concatenate packets
packetsConc = tmri.concatenatePackets(thePackets);
% put HRF back in
packetsConc.HRF = thePackets{1}.HRF;
% locking matrix
paramLockMatrix = createParamLockMatrixVanilla(unique(packetsConc.metaData.theFrequencyIndices) ...
                                               ,packetsConc.metaData.theFrequencyIndices,2);

defaultParamsInfo.nEvents = size(packetsConc.stimulus.values,1);

%% Set parameters
%
% Six parameters define a quadratic form in three dimensions, but
% we normalize the first to 1 so we only need five numbers here.
params0 = tmri.defaultParams('DefaultParamsInfo',defaultParamsInfo);
fprintf('Default model parameters:\n');
tmri.print(params0);

%% Test paramsToVec and vecToParams
params1 = params0;
params1.paramsMainMatrix = rand(size(params1.paramMainMatrix));
x1 = tmri.paramsToVec(params1);
params2 = tmri.vecToParams(x1);
if (any(params1.paramMainMatrix ~= params2.paramMainMatrix))
    error('vecToParams and paramsToVec do not invert');
end
%%
% do the fit
[paramsFit,fVal,fitResponse] = tmri.fitResponse(packetsConc,'DefaultParamsInfo',defaultParamsInfo, ...
                          'paramLockMatrix',paramLockMatrix);
fprintf('Model parameter from fits:\n');

% plot amplitudes
[~, meanParamValues,stdErrorParamValues] = tmri.plotParams(paramsFit,packetsConc.metaData.theFrequencyIndices);
theFreq = thePackets{1}.stimulus.metaData.params.theFrequenciesHz ...
           (2:length(thePackets{1}.stimulus.metaData.params.theFrequenciesHz));
theFreqForPlot = [1 theFreq];
figure;
plot(theFreqForPlot,meanParamValues(1,:),'-kd','LineWidth',2,'MarkerSize',15); hold on
set(gca,'Xscale','log'); xlabel('Temporal Frequency (Hz)'); ylabel('Amplitude');
set(gca,'FontSize',15); set(gca,'Xtick',theFreq); axis square;
% plot fit
tmri.plot(packetsConc.response.timebase,packetsConc.response.values); hold on;
tmri.plot(packetsConc.response.timebase,fitResponse,'Color',[0 1 0],'NewWindow',false);
