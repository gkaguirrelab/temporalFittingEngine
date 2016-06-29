function [HRF, cleanedData] = fitHRF_FIR(timeSeries,attnStartTimes,lengthAttnHRF,T_R)

% function HRF = fitHRF_FIR(timeSeries,attnStartTimes)
%
% derives HRF across multiple runs

% call the HRF generator for each run, then average across runs
for i = 1:size(timeSeries,1)
   [hrf,timeSeriesNoAttn] = attentionFIR(T_R.*(1:length(timeSeries(i,:))), ...
             timeSeries(i,:), attnStartTimes(i,attnStartTimes(i,:)~=-1)', ...
             lengthAttnHRF,T_R) ;
   hrfStore(i,:) = hrf ;    
   cleanedData(i,:) = timeSeriesNoAttn;
end

HRF = mean(hrfStore);

gribble = 1;