function U = arap_positions(V,F,R)
% Given a local rotation for each vertex in the mesh),
% this algorithm calculates new positions that minimize
% the ARAP energy (also known as global step in the ARAP method)
% subject to fixing the positions of the vertices of the triangles 
% with lowest barycenter z-coordinate 
%
% Example usage:
% U = arap_positions(V,F,R)
%
% Input:
% V: #vertices by 3 list of vertex positions
% F: #faces by 3 list of face vertex indices
% R: 3 by 3 by #vertices list of per-vertex given rotations
%
% Output:
% U: #vertices by 3 list of new vertex positions
%

% Bary = barycenter(V,F);
% [~,fixed_tri] = min(Bary(:,3));

% parameters to determine fixed points
cells_per_dim = 2;
xmin = min(V(:,1));
xmax = max(V(:,1));
xlen = (xmax-xmin)/cells_per_dim;
ymin = min(V(:,2));
ymax = max(V(:,2));
ylen = (ymax-ymin)/cells_per_dim;
zmin = min(V(:,3));
zmax = max(V(:,3));
zlen = (zmax-zmin)/cells_per_dim;

% determining fixed vertices
fixed_V = zeros(cells_per_dim,cells_per_dim,cells_per_dim);
for i=1:cells_per_dim
    for j=1:cells_per_dim
        for k=1:cells_per_dim
            verts = find((V(:,1)>(xmin+(i-1)*xlen))&(V(:,1)<(xmin+i*xlen))...
                &(V(:,2)>(ymin+(j-1)*ylen))&(V(:,2)<(ymin+j*ylen))...
                &(V(:,3)>(zmin+(k-1)*zlen))&(V(:,3)<(zmin+k*zlen)));
            
            cell_center = [xmin+(i-1/2)*xlen ymin+(j-1/2)*ylen zmin+(k-1/2)*zlen];
            min_dist = inf;
            for a=1:size(verts,1)
                if (norm(verts(a,:)-cell_center))<min_dist
                    min_dist = norm(V(verts(a),:)-cell_center);
                    fixed_V(i,j,k) = verts(a);
                end
            end
        end
    end
end
fixed_V = setdiff(fixed_V,[0]);
% tsurf(F,V);
% axis equal
% hold on
% for a=1:size(fixed_V,1)
%     plot3(V(fixed_V,1),V(fixed_V,2),V(fixed_V,3),'.r','MarkerSize',20)
% end

L = cotmatrix(V,F);
B = arap_rhs(V,F,R);
zQ = -0.5*L;
zL = -B;

U = min_quad_with_fixed(zQ,zL,fixed_V,V(fixed_V,:));
% obs. 1: decide if it necessary to fix some position
% (as it is now, the result may contain a global translation)
% obs. 2: the unconstrained problem is equivalent to solving
% U = (-L)\B;
