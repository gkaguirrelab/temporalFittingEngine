% Examine parameter fragility
%
% How robust is the partition of variance between the transient and
% sustained components? Is this partition recovered in the face of noise
% and variations in the gamma component?
%


close all
clear all

%lms
%params=[0.1010    0.2213    0.1963    0.8305    0.1618    0.0006    1.3247];

% mel
params=[0.1020    0.2257    0.0915    0.9459    0.0649    0.5063    1.1021];

t=linspace(0,13,13000)';
noise=randn(130,1)*.005;
noise=interp(noise,100);
stimulus=createPupilStimulusVector();

adjustment=linspace(0.75,1.25,10);

% adjust the relative amplitude of the transient and sustained components
for f=1:length(adjustment)
    adjParams=params;
    adjParams(3)=adjParams(3)*adjustment(f);
    adjParams(5)=adjParams(5)/adjustment(f);
    [yPupil,~,~,~] = createPupilTemporalModel(t,stimulus,adjParams);
    [~, fitParams, ~] = fitPupilModelToData(t, yPupil+noise, noise, stimulus);
    result(1,f)=(fitParams(3)/fitParams(5))/(params(3)/params(5));
end

% adjust the time constant of the gamma and see if the sustained amplitude
% stays the same

for f=1:length(adjustment)
    adjParams=params;
    adjParams(2)=adjParams(2)*adjustment(f);
    [yPupil,~,~,~] = createPupilTemporalModel(t,stimulus,adjParams);
    [~, fitParams, ~] = fitPupilModelToData(t, yPupil+noise, noise, stimulus);
    result(2,f)=fitParams(3)/params(3);
end


% Plot the results
hold off
figure
r1=plot(adjustment,result(1,:),'-o');
hold on
r2=plot(adjustment,result(2,:),'-o');
xlabel('Adjustment ratio');
ylabel('Result');
legend([r1 r2], 'Amp sustained:persistent', 'Adjust Gamma, measure amp sustained'); legend boxoff;

