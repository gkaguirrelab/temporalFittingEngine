function paramLockMatrix = daisyChainParamLock(numParam,params2lock)

% function paramLockMatrix = daisyChainParamLock(numParam,params2lock)
%
% daisy chains parameter matrix for locking
%
% numParam   : total number of parameters (number of columns in locking
%              matrix)
% params2lock: matrix with each row specifying a set of parameters
%              (expressed in terms of their positions in the matrix) to be
%              chained together. Since each row of this matrix might have a
%              different length, the ends of each row should be padded by
%              an appropriate number of -1 values

paramLockMatrix = [];

% for each set of parameters to chain together
for i = 1:size(params2lock,1)
    % get those parameters
    params2chain = params2lock(i,:);
    params2chain = params2chain(params2chain>0);
    % chain each successive pair together
   for j = 1:length(params2chain)-1
       % initialize row
       constraintVec = zeros([1 numParam]); 
       % first value in daisy chain pair = 1, second = -1; that is, first
       % parameter - second parameter = 0. i.e., they are equal
       constraintVec(params2chain(j)) = 1;
       constraintVec(params2lock(j+1)) = -1;
       paramLockMatrix(size(paramLockMatrix,1),:) = constraintVec;
   end
end

gribble = 1;