function varargout = uttb_tfeApplyKernel(varargin)
% varargout = uttb_tfeApplyKernel(varargin)
%
% Works by running v_tfeApplyKernel with various arguments and comparing
% results with those stored.
%
% Validate applyKernel method of tfe parent class.

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)
    
    
    %% Basic validation
    UnitTest.validationRecord('SIMPLE_MESSAGE', '***** v_tfeApplyKernel *****');
    validationData1 = v_tfeApplyKernel('makePlots',runTimeParams.generatePlots);
    UnitTest.validationData('validationData1',validationData1);
    
    %% Change kernel timebase deltaT
    UnitTest.validationRecord('SIMPLE_MESSAGE', '***** v_tfeApplyKernel(''kernelDeltaT'',0.5) *****');
    validationData2 = v_tfeApplyKernel('kernelDeltaT',0.5,'makePlots',runTimeParams.generatePlots);
    UnitTest.validationData('validationData2',validationData2,'makePlots',runTimeParams.generatePlots);
    
    %% Change kernel timebase deltaT to an odd value
    UnitTest.validationRecord('SIMPLE_MESSAGE', '***** v_tfeApplyKernel(''kernelDeltaT'',0.8) *****');
    validationData2 = v_tfeApplyKernel('kernelDeltaT',0.8,'makePlots',runTimeParams.generatePlots);
    UnitTest.validationData('validationData2',validationData2,'makePlots',runTimeParams.generatePlots);
    
end



