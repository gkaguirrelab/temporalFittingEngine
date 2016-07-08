
%% Fake Data to test Fourier & FIR

Time_Samples = 1:336 ;
Fake_Data = zeros(1,336) ;

AttnStartTimes = [% 21.6630000000000;...
%29.5720000000000;77.9870000000000;139.510000000000;...
159.030000000000];%224.090000000000;234.130000000000;...
%291.280000000000;308.820000000000] ; %320.610000000000;...
%333.530000000000] ;

Rndnumbers = [%50; 90; ..
    160]; %230; 270; 300] ;

AttnRound = floor(AttnStartTimes) ;

for i = 1:length(Rndnumbers)
    Temp = Rndnumbers(i) ;
    Fake_Data(Temp:Temp + 10) = 1 ;
end 

BOLDHRF = createCanonicalHRF(Time_Samples,6,12,10);
% interp1(Fake_Data,BOLDHRF) ;

Fake_TimeSeries = conv(Fake_Data,BOLDHRF) ;
Fake_TimeSeries = Fake_TimeSeries(1:336) ;



[hrf] = attentionFourier(Time_Samples,Fake_TimeSeries,AttnStartTimes,16,1) ;
hrf_Fourier = hrf ;
clear hrf ;

[hrf] = attentionFIR(Time_Samples,Fake_TimeSeries,AttnStartTimes,16,1) ;
hrf_FIR = hrf ;

figure ;
subplot(2,1,1)
plot(hrf_FIR) ;
axis square ;
subplot(2,1,2)
plot(hrf_Fourier) ;
axis square ;


Graves = 69 ;