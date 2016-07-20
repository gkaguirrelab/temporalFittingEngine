function tmriPrintParams(obj,varargin)
% tmriPrintParams(obj,varargin)
%
% Print out the parameters to the command window
%
% Keys
%   'ParamType',['simulate'/'fit'] - which parameters to print 
%
% Can print either simulate or fit parameters

%% Parse string value pairs
% p = inputParser;
% p.addParameter('ParamType','simulate',@isstring);
% p.parse(varargin{:});

% switch (p.Results.ParamType)
%     case 'simulate'
%         params = obj.simulateParams;
%     case 'fit'
%         params = obj.fitParams;
%     otherwise
%         error('Unknown parameter type specified');
% end
params = obj.simulateParams;

% Quadratic parameters
fprintf('Quadratic ellipse lengths: 1, %0.2f, %0.2f\n',params.Qvec(1),params.Qvec(2));
fprintf('Quadratic ellipse angles (degs): %0.1f, %0.1f %0.1f\n',params.Qvec(3),params.Qvec(4),params.Qvec(5));
fprintf('CRF amplitude: %0.2f, CRF semi-saturation: %0.2f, CRF exponent: %0.2f\n',params.crfAmp,params.crfSemi,params.crfExponent);
fprintf('Exponential filter time constant: %0.2f\n',params.expFalloff);

end