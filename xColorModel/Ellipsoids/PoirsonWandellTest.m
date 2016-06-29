% PoirsonWandellTest
%
% Test the routines that generate data according to the Poirson-Wandell
% model.
%
% 6/29/16  dhb  Wrote it.

%% Clear and close
clear; close all;

conditionStr = 'HT,cc';
theSf = 2;
[A,Ainv,Q] = PoirsonWandellEllipsoid(conditionStr,theSf);

nTheta = 20;
nPhi = 20;
xSphere = GeneratePointsOn3DUnitSphere(nTheta,nPhi);
xEllipsoid = Ainv*xSphere;

% Plot the fit as a nice surface
xCoords = squeeze(xEllipsoid(1,:));
yCoords = squeeze(xEllipsoid(2,:));
zCoords = squeeze(xEllipsoid(3,:));
tri = delaunay(xCoords, yCoords, zCoords);

figure; clf;
subplot(1,3,1); hold on
h = trisurf(tri, xCoords, yCoords, zCoords);
set(h,'FaceAlpha',0.25)
set(h,'EdgeColor',[0.5 0.5 0.5])
set(h,'FaceColor',[0.6 0.6 0.6]);
lighting phong;
xlabel('L contrast'); ylabel('M contrast'); zlabel('S contrast'); title('Poirson Wandell Ellipsoid');
%xlim([-0.15 0.15]); ylim([-0.15 0.15]); zlim([-0.25 0.25]);
axis('square');

subplot(1,3,2); hold on
h = trisurf(tri, xCoords, yCoords, zCoords);
set(h,'FaceAlpha',0.25)
set(h,'EdgeColor',[0.5 0.5 0.5])
set(h,'FaceColor',[0.6 0.6 0.6]);
lighting phong;
xlabel('L contrast'); ylabel('M contrast'); zlabel('S contrast'); title('Poirson Wandell Ellipsoid');
%xlim([-0.15 0.15]); ylim([-0.15 0.15]); zlim([-0.25 0.25]);
axis('square');
view(0,0);

subplot(1,3,3); hold on
h = trisurf(tri, xCoords, yCoords, zCoords);
set(h,'FaceAlpha',0.25)
set(h,'EdgeColor',[0.5 0.5 0.5])
set(h,'FaceColor',[0.6 0.6 0.6]);
lighting phong;
xlabel('L contrast'); ylabel('M contrast'); zlabel('S contrast'); title('Poirson Wandell Ellipsoid');
%xlim([-0.15 0.15]); ylim([-0.15 0.15]); zlim([-0.25 0.25]);
axis('square');
view(90,0);