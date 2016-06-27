function xEllipsoid = GenerateEllipsoid(ellParams,nTheta,nPhi)
% xEllipsoid = GenerateEllipsoid(ellParams,nTheta,nPhi)
%
% Generate a set of points on the 3D ellipsoid by stretching and rotating
% points on the 3D unit sphere.
%
% See GenerateEllipsoidMatrices for specification of ellParams vector.
%
% 6/27/16  dhb  Wrote it.

xSphere = GeneratePointsOn3DUnitSphere(nTheta,nPhi);
[A,Ainv,Q] = GenerateEllipsoidMatrices(ellParams);
xEllipsoid = Ainv*xSphere;
