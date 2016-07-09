
%% Fake Data to test Fourier & FIR

Time_Samples = 1:336 ;
Fake_Data = zeros(1,336) ;

AttnStartTimes = [ 21.6630000000000;...
29.5720000000000;77.9870000000000;139.510000000000;...
159.030000000000; 224.090000000000;234.130000000000;...
291.280000000000;308.820000000000; 320.610000000000;...
333.530000000000] ;

% AttnStartTimes = [224.09];

AttnRound = round(AttnStartTimes) ;

Rndnumbers = AttnRound;

for i = 1:length(Rndnumbers)
    Temp = Rndnumbers(i) ;
    Fake_Data(Temp) = 1 ;
end 

gamma1 = 6; gamma2 = 12; gammaScale = 10; HRFexpectedLength = 16;

BOLDHRF = createCanonicalHRF(Time_Samples,gamma1,gamma2,gammaScale);
BOLDHRFtoTest = createCanonicalHRF(1:HRFexpectedLength,gamma1,gamma2,gammaScale);
% interp1(Fake_Data,BOLDHRF) ;

Fake_TimeSeries = conv(Fake_Data,BOLDHRF) ;
Fake_TimeSeries = Fake_TimeSeries(1:336) ;

[hrf] = attentionFourier(Time_Samples,Fake_TimeSeries,AttnStartTimes,HRFexpectedLength,1) ;
hrf_Fourier = hrf ;
hrf_Fourier = hrf_Fourier - hrf_Fourier(1);
clear hrf ;

[hrf] = attentionFIR(Time_Samples,Fake_TimeSeries,AttnStartTimes,HRFexpectedLength,1) ;
hrf_FIR = hrf ;
hrf_FIR = hrf_FIR - hrf_FIR(1);

figure ;
set(gcf,'Position',[401 418 1474 460]);

subplot(2,2,2)
plot(hrf_FIR) ; title('FIR'); axis square;

subplot(2,2,1)
plot(1:HRFexpectedLength,BOLDHRFtoTest); axis square;
title('Ground-truth HRF');

subplot(2,2,3)
plot(Time_Samples,Fake_TimeSeries); title('Fake Timeseries');

subplot(2,2,4)
plot(hrf_Fourier) ; title('Fourier');
axis square ;


Graves = 69 ;