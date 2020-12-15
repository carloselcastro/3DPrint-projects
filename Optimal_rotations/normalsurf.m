function N = normalsurf(V,F)
%Find the Normal field of surface V and triangle index F.
%
%Sintax:
%
%   N = normalsurf(V,F)
%
% Inputs:
%
%   V       #V x 3 matrix of surface's vertex coordinates.
%   F       #F x 3  matrix of indices of surface's triangle corners.

% Output:
%   N       #F x 3 matrix of Normal field of surface.
%

v = V(F(:,3),:) - V(F(:,1),:); % Setting triangle's edges
w = V(F(:,2),:) - V(F(:,1),:); % Setting triangle's edges
vxw=cross(v,w); % Cross product of triangle's edges and finding normals
N=normr(vxw); % Normalizing the normals
end