fullFigure;
numFreqs = [16 25 33 35];
scaleFactor = 10000;
for i = 1:length(numFreqs)
    sampT                   = 1000; % TR (msec)
    eventTimes              = [12200 17325 33300 55500 60320 88800]; % attention events (msec)
    HRFdur                  = numFreqs(i) * 1000; % msec 
    %HRFdur                  = 30000; % msec
    % Make fake data
    hrf                     = doubleGammaHrf(1/sampT,[6 10],[1 1],0,HRFdur/1000); % HRF in sec
    hrf                     = hrf * scaleFactor;
    tmpTC                   = zeros(100000,1);
    tmpTC(eventTimes)       = 1;
    upTC         = filter(hrf,1,tmpTC);
    % add some noise
    nupTC                   = upTC + 0*max(upTC)*randn(size(upTC));
    timeSeries              = downsample(nupTC,sampT);
    %
    subplot(1,4,i);
    [HRF, timeSeriesNoAttn] = attentionFIR(0:length(timeSeries)-1,timeSeries'-mean(timeSeries),round(eventTimes./sampT)',HRFdur./sampT,1);
    % Plot
    x = (0:length(HRF)-1).*1000;
    plot(x,HRF,'b',0:length(hrf)-1,hrf,'r');
    xlim([0 35000]);
    axis square;
    xlabel('Time (TRs)','FontSize',20);
    ylabel('Percent Signal Change','FontSize',20);
    legend({'FIR model' 'Simulated data'},'FontSize',10);
%    title(['nFreqs = ' num2str(numCov) ],'FontSize',20);
end

%%
%% Set params
sampT                       = 1000; % TR (msec)
numFreqs                    = [32];
scaleFactor                 = 10000; % scaling factor (for HRF
numEvents                   = 10;
offTR = 1; % offTR = 1 puts the events in the middle of the TR
overlap = 5; % overlap = N duplicates the events, spaced apart by N TRs
%%
fullFigure;
for i = 1:length(numFreqs)
    tcDur                   = sampT * numFreqs(i) * 100; % total time-series duration
    trTimes                 = 0:sampT:tcDur-sampT; % TR start times
    eventSpace              = size(trTimes,2)/numEvents; % event spacing
    eventTimes              = trTimes( (eventSpace:eventSpace:size(trTimes,2)) - (eventSpace-1));
    if offTR
        eventTimes = eventTimes + sampT/2;
    end
    if overlap
        eventTimes = sort([eventTimes,eventTimes + overlap*sampT]);
    end
    HRFdur                  = numFreqs(i) * 1000; % msec
    % Make fake data
    hrf                     = doubleGammaHrf(1/sampT,[6 10],[1 1],1/6,HRFdur/1000); % HRF in sec
    hrf                     = hrf * scaleFactor;
    % Make hrf using a sinusoid
    %     t                       = linspace(0,HRFdur-1,HRFdur);
    %     hrf                     = sin(t/HRFdur*2*pi);
    tmpTC                   = zeros(tcDur,1);
    tmpTC(eventTimes+1)     = 1; % add one for indexing
    upTC                    = filter(hrf,1,tmpTC);
    % add some noise
    %nupTC                   = upTC + 0.05*max(upTC)*randn(size(upTC));
    timeSeries              = resample(upTC,1,sampT);
    % Plot upsampled HRFs
    subplot(2,length(numFreqs),i);
    [HRF,fSet,betaValues,DesignMatrix,numCov] = deriveHRF(timeSeries,eventTimes,sampT,HRFdur,numFreqs(i));
    disp(numCov);
    x = 1:length(HRF);
    plot(x,HRF,'b',x,hrf,'r');
    xlim([0 HRFdur]);
    axis square;
    xlabel('Time (msec)','FontSize',20);
    ylabel('Percent Signal Change','FontSize',20);
    legend({'Fourier model' 'Simulated data'},'FontSize',10);
    title(['HRF duration = ' num2str(HRFdur/1000) 's'],'FontSize',20);
    % Plot simulated data
    subplot(2,length(numFreqs),i+length(numFreqs));
    plot(timeSeries);
    xlim([0 length(timeSeries)]);
    axis square;
    xlabel('Time (seconds)','FontSize',20);
    ylabel('Percent Signal Change','FontSize',20);
    title('Simulated data','FontSize',20);
    % Plot downsampled HRFs
    %     subplot(3,length(numFreqs),i+2*length(numFreqs));
    %     dx = 1:max(x)/sampT;
    %     dHRF = resample(HRF,1,sampT);
    %     dhrf = resample(hrf,1,sampT);
    %     plot(dx,dHRF,'b',dx,dhrf,'r');
    %     xlim([0 max(dx)]);
    %     axis square;
    %     xlabel('Time (sec)','FontSize',20);
    %     ylabel('Percent Signal Change','FontSize',20);
    %     legend({'Fourier model' 'Simulated data'},'FontSize',10);
    %     title(['HRF duration = ' num2str(HRFdur/1000) 's'],'FontSize',20);
end
%%
for i = 1:length(betaValues)
    tmp(i,:) = betaValues(i) * fSet(:,i);
end
%%
figure;
tmpSet = fSet(:,2:end);
for i = 1:size(tmpSet,2);
    subplot(7,10,i);
    plot(tmpSet(:,i));
    xlim([0 size(tmpSet,1)]);
end