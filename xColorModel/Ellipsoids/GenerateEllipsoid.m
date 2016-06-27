function x = GenerateEllipsoid(nTheta,nPhi,mean,Q)
% x = GenerateEllipsoid(nPoints,mean,Q)
%
% Generate a set of points on the 3D ellipsoid defined as follows.
%   1) Generate a set of (nTheta*nPhi) points m on the 3D unit sphere using
%   GeneratePointsOn3DUnitSphere.
%   2) Find the mapping between m and x, such that x'*Q*x = 1.
%     a) Write Q = A'*A
%     b) Then x = A'*m

