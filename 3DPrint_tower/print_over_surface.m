function [X,Vnew,Fnew,minim] = print_over_surface(V,F)
%    This algorithm searches for the best orientation in
% space that solves the minimization problem:
%
%                   min sum W*||N-PrN||, 
%                   
% where the weight W is determined using the height to a copy 
% that was printed below.
%
%Syntax:
%
% [X,Vnew,F,fval,exitflag,output] = print_over_surface(V,F)
%
%Input:
%
%   V           #V by 3 vertex positions of the mesh to be printed
%   F           #F by 3 vertex indices of each triangle of the mesh to be
%               printed
%
%Outputs:
%
%   Xmin        2 by 1 vector solution of the minimization problem
%   Vnew        #V by 3 list of vertices at the optimal position
%   F           #F by 3 list of triangle indices
%   minim       objective function value at the optimal point


% 1) Preprocessing

% translate the first copy to a position below the second copy such that
% the two copies don't intersect: (V_prev, F_prev) is going to be the copy 
% printed below and (V,F) is going to be the copy to be printed above
z_min = min(V(:,3));
z_max = max(V(:,3));
V_prev(:,1:2) = V(:,1:2);
V_prev(:,3) = V(:,3)-(z_max-z_min);
F_prev = F;
% translate both so that (V,F) is centred at the origin (to perform
% rotations around the origin)
Xbari = (V(F(:,1),1) + V(F(:,2),1) + V(F(:,3),1))/3;
Ybari = (V(F(:,1),2) + V(F(:,2),2) + V(F(:,3),2))/3;
Zbari = (V(F(:,1),3) + V(F(:,2),3) + V(F(:,3),3))/3;
X = sum(Xbari)/length(Xbari);                                 
Y = sum(Ybari)/length(Ybari);                                 
Z = sum(Zbari)/length(Zbari);
V(:,1) = V(:,1)-X;
V(:,2) = V(:,2)-Y;
V(:,3) = V(:,3)-Z;
V_prev(:,1) = V_prev(:,1)-X;
V_prev(:,2) = V_prev(:,2)-Y;
V_prev(:,3) = V_prev(:,3)-Z;
% visualization
h_prev = tsurf(F_prev,V_prev); axis equal;
hold on
tsurf(F,V); axis equal;
input('')

% 2) sample height field defined by (V_prev,F_prev)
nx = 50;
ny = 50;
sampled_HF = sample_height_field(V_prev,F_prev,nx,ny);
h_points = plot3(sampled_HF(:,:,1),sampled_HF(:,:,2),sampled_HF(:,:,3),'.r','markerSize',10)
input('')
% visualization
delete(h_prev);
delete(h_points);
surf(sampled_HF(:,:,1),sampled_HF(:,:,2),sampled_HF(:,:,3),'FaceAlpha',0.2)

% 3) (DONE) Write auxiliary functions interpolate_height_field and 
% vertical_distance_to_height_field

% 4) print3Dopt_grid with extra input: ('zmin', sampled_HF)).
%
% Replace
%
% (-sin(x(2,1))*...                                               %
%         Bary(:,1)+ cos(x(2,1))*sin(x(1,1))*Bary(:,2) + ...         %
%           cos(x(1,1))*cos(x(2,1))*Bary(:,3) - min(-sin(x(2,1))...  %
%           *Bary(:,1)+ cos(x(2,1))*sin(x(1,1))*Bary(:,2) + ...      %
%           cos(x(1,1))*cos(x(2,1))*Bary(:,3)))                      %
% 
% which are the z coordinates of the rotated barycenters minus z_min 
%
% by the vertical distance of the rotated barycenters to the sampled 
% height field. The expressions for the 3 coordinates of the rotated 
% barycenters will have to be explicitly calculated. Then use 
%
% vertical_distance_to_height_field("expression for rotated barycenters", sampled_HF)


% placeholder outputs
X = [];
Vnew = [];
Fnew = []; 
minim = [];