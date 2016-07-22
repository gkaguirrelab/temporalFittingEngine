classdef tmriQuadraticModelTest < matlab.unittest.TestCase
    
    properties
        testFigure;
    end
    
    methods (TestClassSetup)
        function initTests(obj)
            AddToMatlabPathDynamically(fullfile(fileparts(which(mfilename)),'..','toolbox'));
            obj.testFigure = figure; clf; hold on
        end
    end
    
    methods (TestClassTeardown)
        function closeTests(obj)
            close(obj.testFigure);
        end
    end
    
    methods (Test)
        % Test that instantiating the class doesn't crash
        function tmriModelConstructorTest(obj)
            tmri = tmriQuadraticColorModel;
        end
        
        % Test that paramsToVec and vecToParams invert
        % 
        % This also tests that defaultParams returns a parameter struct and
        % that print prints it out.
        function paramsToVecTest(obj)
            tmri = tmriQuadraticColorModel;
            params0 = tmri.defaultParams;
            tmri.print(params0);
            
            x0 = tmri.paramsToVec(params0);
            x1 = x0;
            x1(1) = 2;
            x1(2) = 0.5;
            x1(3) = pi/2;
            x1(7) = 3;
            params1 = tmri.vecToParams(x1);
            x2 = tmri.paramsToVec(params1);
            obj.assertEqual(x1,x2);
        end
        
        % Test that we can simulate a neural response
        function neuralResponseTest(obj)
            % Construct the model object
            tmri = tmriQuadraticColorModel;
            
            % Set parameters
            params0 = tmri.defaultParams;

            % Set the timebase we want to compute on
            deltaT = 1;
            totalTime = 1000;
            timebase = 0:deltaT:totalTime;
            
            % Specify the stimulus.
            nTimeSamples = size(timebase,2);
            filter = fspecial('gaussian',[1 nTimeSamples],6);
            stimulus= rand(3,nTimeSamples);
            for i = 1:3
                stimulus(i,:) = ifft(fft(stimulus(i,:)) .* fft(filter));
            end
            
            % Generate response and plot
            neuralResponse = tmri.computeResponse(params0,timebase,stimulus,'AddNoise',true);
            figure(obj.testFigure);
            tmri.plot(timebase,neuralResponse,'NewWindow',false);
        end
    end
    
end