% EllipsoidTest
%
% Test the ellipsoid code
%
% 6/27/16  dhb  Wrote it.

%% Clear
clear; close all;

%% Generate points on unit sphere
nTheta = 20;
nPhi = 20;
xSphere = GeneratePointsOn3DUnitSphere(nTheta,nPhi);
figure; clf; hold on
plot3(xSphere(1,:),xSphere(2,:),xSphere(3,:),'ro','MarkerSize',8,'MarkerFaceColor','r');
axis('square');
xlabel('X'); ylabel('Y'); zlabel('Z'); title('Unit Sphere');

%% Generate the ellipsoid matrices from a vector
% of ellipsoid parameters
ellParams = [4 0.5 2 pi/4 pi/8 pi]';
[A,Ainv,Q] = GenerateEllipsoidMatrices(ellParams);

%% Map the unit sphere into the ellipsoid
xEllipsoid = Ainv*xSphere;
figure; clf; hold on
plot3(xEllipsoid(1,:),xEllipsoid(2,:),xEllipsoid(3,:),'ro','MarkerSize',8,'MarkerFaceColor','r');
xlabel('X'); ylabel('Y'); zlabel('Z'); title('Ellipsoid');

%% Make sure individual ellipsoid points satisfy xEllipsoid'*Q*xEllipsoid = 1
check = diag(xEllipsoid'*Q*xEllipsoid);
tolerance = 1e-8;
if (max(abs(check(:)-1)) > tolerance)
    error('Failure of ellipsoid points to have transformed length of 1');
end

%% Map ellipsoid back into unit sphere space and make sure it is a sphere
xSphereCheck = A*xEllipsoid;
figure; clf; hold on
plot3(xSphereCheck(1,:),xSphereCheck(2,:),xSphereCheck(3,:),'ro','MarkerSize',8,'MarkerFaceColor','r');
axis('square');
xlabel('X'); ylabel('Y'); zlabel('Z'); title('Unit Sphere Check');

