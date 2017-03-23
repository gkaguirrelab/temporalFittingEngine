function varargout = v_IAMPBasic(varargin)
% function varargout = v_IAMPBasic(varargin)
%
% Works by running t_IAMPBasic with various arguments and comparing
% results with those stored.
%

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)
    
    %% Freeze RNG so validations work
    rng(1);
    
    %% Basic validation
    UnitTest.validationRecord('SIMPLE_MESSAGE', '***** v_IAMPBasic *****');
    validationData1 = t_IAMPBasic('generatePlots',runTimeParams.generatePlots);
    UnitTest.validationData('validationData1',validationData1, ...
        'UsingTheFollowingVariableTolerancePairs', ...
         'validationData1.paramMainMatrix',5e-4);
    
end
