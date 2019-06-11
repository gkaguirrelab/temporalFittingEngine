function joinedPacket = joinPackets(obj,packetsCellArray)
% Joins the subfields of packets together 
%
% Syntax:
%    joinedPacket = OBJ.joinPackets(packetsCellArray)
%
% Description:
%    This function takes in a cell array of packets an join the subfields 
%    together to create one packet where each subfield is N times longer.
%    This keeps the same number of  regressors for the IAMP model when the 
%    pakets are joined (this is not what the current concatenatePacket method
%    does in tfe, it creates N x numRegressors). 
%
% Inputs:
%    packetsCellArray          - Cell array of packets to be joined
%
% Outputs:
%    joinedPacket              - The joined packets 
% Optional key/value pairs:
%    none

% MAB 06/11/19

%% Initialize dummy packet
joinedPacket.response.values   = [];
joinedPacket.response.timebase = [];
% the stimulus
joinedPacket.stimulus.timebase = [];
joinedPacket.stimulus.values   = [];
% the kernel
joinedPacket.kernel = [];
% the meta data (this is the constrast and directions)
joinedPacket.metaData = [];


for ii = 1:length(packetsCellArray)
    
    % the response
    joinedPacket.response.values   = [joinedPacket.response.values packetsCellArray{ii}.response.values];
    
    % the stimulus
    joinedPacket.stimulus.values   = [joinedPacket.stimulus.values packetsCellArray{ii}.stimulus.values];
  
    % Get time base info for all packets
    globalTimebase(ii,:) = packetsCellArray{ii}.response.timebase;
end

deltaT = median(diff(globalTimebase'));
if sum(sum(diff(globalTimebase') ~= median(deltaT))) > 0
    error('Inconsistent timebase')
else
    deltaT = median(deltaT);
end

concatTimebase= deltaT:deltaT:deltaT*length(joinedPacket.response.values);

joinedPacket.response.timebase = concatTimebase;
joinedPacket.stimulus.timebase = concatTimebase;
% the kernel
joinedPacket.kernel = generateHRFKernel(6,12,10,concatTimebase);
% the meta data (this is the constrast and directions)
joinedPacket.metaData = [];

end


