function out = padResample(X,P,Q,padType)

%   Resamples the values, X, of a uniformly sampled signal at P/Q times
%   the original sample rate using a polyphase anti-aliasing filter.
%   This is doen by calling the 'resample' function, but first padding the
%   start and end using data taken from X. The type of padding is specified
%   by the 'padType' input.
%
%   Motivation for this function comes from the 'resample' function header:
%
%   In its filtering process, RESAMPLE assumes the samples at times before
%   and after the given samples in X are equal to zero. Thus large
%   deviations from zero at the end points of the sequence X can cause
%   inaccuracies in Y at its end points.
%
%   Usage:
%   out = padResample(X,P,Q,padType)
%
%   Inputs:
%   X       = input data matrix
%   P       = numerator of output sample rate
%   Q       = denominator of output sample rate
%   padType = type of padding (default = 'periodic')
%
%   padType options:
%   'periodic'      = [X;X;X]; <default>
%   'mirror'        = [flipud(X);X;flipud(X)];
%   'firstLast'     = [repmat(X(1,:),size(X,1),1);X;repmat(X(end,:),size(X,1),1)];
%   'zero'          = [zeros(size(X));X;zeros(size(X))];
%
%   Note:
%   Assumes the matrix X is organized such that size(X,1) > size(X,2)
%
%   Written by Andrew S Bock Aug 2016

%% set defaults
if ~exist('padType','var')
    padType = 'periodic';
end
%% Check the size of X
if size(X,2) > size(X,1)
    error('size(X,2) > size(X,1), consider transposing X');
end
%% Pad the input using first and last values
% padStart    = ones(size(X)) .* repmat(X(1,:),size(X,1),1);
% padEnd      = ones(size(X)) .* repmat(X(end,:),size(X,1),1);
% padX        = [padStart;X;padEnd];
%% Pad the input using fliplr(X) and X
switch padType
    case 'periodic'
        padX        = [X;X;X];
    case 'mirror'
        padX        = [flipud(X);X;flipud(X)];
    case 'firstLast'
        padX        = [repmat(X(1,:),size(X,1),1);X;repmat(X(end,:),size(X,1),1)];
    case 'zero'
        padX        = [zeros(size(X));X;zeros(size(X))];
end
%% Call resample
padXr       = resample(padX,P,Q);
%% Output the central, non-padded portion
out = padXr(size(X,1)/Q+1:end-size(X,1)/Q,:);