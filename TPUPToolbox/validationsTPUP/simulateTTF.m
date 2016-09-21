
% Simulate TTF
%
% Using the fit parameters for subject HERO_gka1, submit sinusoidal
% modulations to the createPupilTemporalModel routine, and record the
% resulting amplitude. Plot the TTF and compare this to the Spitschan 2014
% PNAS paper result.

close all
clear all

lmsParams=[0.1010    0.2213    0.1963    0.8305    0.1618    0.0006    1.3247];
melParams=[0.1020    0.2257    0.0915    0.9459    0.0649    0.5063    1.1021];

freqs=[0.01,0.05,0.1,0.5,1,2];
t=linspace(0,119,12000)';


% LMS
figure
hold on
for f=1:length(freqs)
    stimulus=sin( (t.*freqs(f))*2*pi );
    [yPupil,~,~,~] = createPupilTemporalModel(t,stimulus,lmsParams);
    plot(t,yPupil);
    lmsAmp(f)=max(yPupil(2000:12000))-min(yPupil(2000:12000));
end
hold off

% Mel
figure
hold on
for f=1:length(freqs)
    stimulus=sin( (t.*freqs(f))*2*pi );
    [yPupil,~,~,~] = createPupilTemporalModel(t,stimulus,melParams);
    plot(t,yPupil);
    melAmp(f)=max(yPupil(2000:12000))-min(yPupil(2000:12000));
end
hold off

figure
semilogx(freqs,lmsAmp,'-k')
hold on
semilogx(freqs,lmsAmp,'xk','MarkerSize',10)
semilogx(freqs,melAmp,'-b')
semilogx(freqs,melAmp,'xb','MarkerSize',10)
xlabel('Freq [Hz]');
ylabel('Amplitude response [proportion change]');
hold off