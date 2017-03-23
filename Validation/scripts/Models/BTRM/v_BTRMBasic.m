function varargout = v_BTRMBasic(varargin)
% function varargout = v_BTRMBasic(varargin)
%
% Works by running t_BTRMBasic with various arguments and comparing
% results with those stored.
%

varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

%% Freeze RNG so validations work
rng(1);

%% Basic validation
UnitTest.validationRecord('SIMPLE_MESSAGE', '***** v_BTRMBasic *****');
validationData = t_BTRMBasic('generatePlots',runTimeParams.generatePlots);
UnitTest.validationData('validationData',validationData, ...
    'UsingTheFollowingVariableTolerancePairs', ...
    'validationData.demo1.paramMainMatrix',5e-4,...
    'validationData.demo3.paramMainMatrix',5e-4,...
    'validationData.demo4.paramMainMatrix',5e-4);

end
