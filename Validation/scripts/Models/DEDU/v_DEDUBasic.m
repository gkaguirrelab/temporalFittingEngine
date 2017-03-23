function varargout = v_DEDUBasic(varargin)
% function varargout = v_DEDUBasic(varargin)
%
% Works by running t_DEDUBasic with various arguments and comparing
% results with those stored.
%

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)
    
    %% Freeze RNG so validations work
    rng(1);
    
    %% Basic validation
    UnitTest.validationRecord('SIMPLE_MESSAGE', '***** v_DEDUBasic *****');
    validationData1 = t_DEDUBasic('generatePlots',runTimeParams.generatePlots);
    UnitTest.validationData('validationData1',validationData1, ...
        'UsingTheFollowingVariableTolerancePairs', ...
         'validationData1.paramMainMatrix',5e-4);
    
end
