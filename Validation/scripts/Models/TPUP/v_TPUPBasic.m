function varargout = v_TPUPBasic(varargin)
% function varargout = v_TPUPBasic(varargin)
%
% Works by running t_TPUPBasic with various arguments and comparing
% results with those stored.
%

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)
    
    %% Freeze RNG so validations work
    rng(1);
    
    %% Basic validation
    UnitTest.validationRecord('SIMPLE_MESSAGE', '***** v_TPUPBasic *****');
    validationData1 = t_TPUPBasic('generatePlots',runTimeParams.generatePlots);
    UnitTest.validationData('validationData1',validationData1, ...
        'UsingTheFollowingVariableTolerancePairs', ...
         'validationData1.params1.paramMainMatrix',5e-4,...
         'validationData1.modelResponseStruct.values',5e-4);

end
