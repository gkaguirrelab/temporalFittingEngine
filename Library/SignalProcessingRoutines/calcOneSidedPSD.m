function [ psdStruct ] = calcOneSidedPSD( dataStruct, varargin )
% function [ psdStruct ] = calcOneSidedPSD( dataStruct )
%
% This routine returns the power spectral denisty of the input signal.
%
% Inputs:
%   dataStruct - a structure with a timebase and values field. If the
%     values field contains a matrix, then a psd will be calculated for
%     each row.
%
% Optional arguments:
%   meanCenter - Logical. Defines if the DC component should be discarded
%   verbosity -
%     'none' - the defaut. Shh.
%     'full'
%
% Outputs:
%   psdStruct - the one-sided power spectrum of the dataStruct. The time-
%     base is set to the one-sided frequency range (in Hz). The length
%     is one-half the input lenght. The values are in units of power, with
%     each row of the values field corresponding to each row of the values
%     field in the input dataStruct. The sum of the values in the one-sided
%     spectrum is equal to the variance of the input signal.
%

%% Parse vargin for options passed here
%
% Setting 'KeepUmatched' to true means that we can pass the varargin{:})
% along from a calling routine without an error here, if the key/value
% pairs recognized by the calling routine are not needed here.
p = inputParser; p.KeepUnmatched = true;
p.addRequired('dataStruct',@isstruct);
p.addParameter('verbosity','none',@ischar);
p.addParameter('meanCenter',false,@islogical);
p.parse(dataStruct, varargin{:});

dataStruct=p.Results.dataStruct;

% check that the dataStruct is well formed
if ~isfield(dataStruct,'values') || ~isfield(dataStruct,'timebase')
    error('The dataStruct does not have both a values and timebase field');
end

dataLength = length(dataStruct.timebase);
nRows=size(dataStruct.values,1);

% check that the length of the timebase is equal to the column length of
% the values field
if dataLength~=size(dataStruct.values,2)
    error('The timebase and values are not the same length');
end

% apologize for not having the solution for odd-length vectors yet
if mod(dataLength,2)
    error('Currently implemented for even-length signals only. Sorry.');
end

% derive the deltaT from the stimulusTimebase (units of msecs)
check = diff(dataStruct.timebase);
deltaT = check(1);

% meanCenter if requested

if p.Results.meanCenter
    dataStruct.values=dataStruct.values - ...
      repmat(mean(dataStruct.values,2),1,dataLength);
end

% Calculate the FFT for each row of the values field.
for ii=1:nRows
    X=fft(dataStruct.values(ii,:));
    psd=X.*conj(X)/(dataLength^2);
    psdStruct.values(ii,:)=psd(1:dataLength/2);
end

% Produce the freq "timebase" in Hz
psdStruct.timebase=deltaT*(0:dataLength/2-1)/dataLength/1000;

% Make a plot if verbose is on
if strcmp(p.Results.verbosity,'full')
    figure

    % plot the time-domain signal
    subplot(2,1,1);
    plot(dataStruct.timebase/1000,dataStruct.values);
    title('Time domain signal');
    xlabel('Time (secs)')
    ylabel('Amplitude');

    % plot the frequency domain signal
    subplot(2,1,2);
    plot(psdStruct.timebase,psdStruct.values,'b','LineWidth',1);
    hold on
    plot(psdStruct.timebase,psdStruct.values,'ro');
    hold off
    title('One Sided Power Spectral Density');
    xlabel('Frequency (Hz)')
    ylabel('Power');
end

end % function