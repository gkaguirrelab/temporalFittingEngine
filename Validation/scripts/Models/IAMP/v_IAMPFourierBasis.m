function varargout = v_IAMPFourierBasis(varargin)
% function varargout = v_IAMPFourierBasis(varargin)
%
% Works by running v_IAMPFourierBasis with various arguments and comparing
% results with those stored.
%

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)
    
    %% Freeze RNG so validations work
    rng(1);
    
    %% Basic validation
    UnitTest.validationRecord('SIMPLE_MESSAGE', '***** v_IAMPFourierBasis *****');
    validationData1 = t_IAMPFourierBasis('generatePlots',runTimeParams.generatePlots);
    UnitTest.validationData('validationData1',validationData1);
    
end
