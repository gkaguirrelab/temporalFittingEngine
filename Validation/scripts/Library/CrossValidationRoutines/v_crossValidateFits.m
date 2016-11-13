function varargout = v_crossValidateFits(varargin)
% function varargout = v_crossValidateFits(varargin)
%
% Works by running t_crossValidateFits with various arguments and comparing
% results with those stored.
%

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)
    
    %% Freeze RNG so validations work
    rng(1);
    
    %% Basic validation
    UnitTest.validationRecord('SIMPLE_MESSAGE', '***** v_crossValidateFits *****');
    validationData1 = t_crossValidateFits('generatePlots',runTimeParams.generatePlots);
    UnitTest.validationData('validationData1',validationData1, ...
        'UsingTheFollowingVariableTolerancePairs', ...
         'validationData1.trainfVals',1e-5, ...
         'validationData1.testfVals',1e-5, ...
         'validationData1.paramMainMatrix',5e-4);
    
end
