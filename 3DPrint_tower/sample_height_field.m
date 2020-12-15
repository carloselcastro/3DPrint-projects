function sampled_HF = sample_height_field(V,F,nx,ny)
%  This functions samples the height field of a given surface (V,F).
%
%Syntax:
%
% sampled_HF = sample_height_field(V,F,nx,ny)
%
%Input:
%
%   V           #V by 3 vertex positions of the mesh to have the height
%               field sampled
%   F           #F by 3 vertex indices of each triangle of the mesh
%               to have the height field sampled
%   nx          number of samples on the x-axis
%   ny          number of samples on the y-axis
%
%Outputs:
%
%   sampled_HF  nx by ny by 3 sampled height field where (:,:,1) stands for
%               x values, (:,:,2) stands for y values and (:,:,3) stands
%               for z (height) values

min_x = min(V(:,1));
max_x = max(V(:,1));
min_y = min(V(:,2));
max_y = max(V(:,2));

min_z = min(V(:,3));

[X,Y] = meshgrid(linspace(min_x,max_x,nx),linspace(min_y,max_y,ny));

% query points
P = zeros(ny*nx,2);
P(:,1) = reshape(X,ny*nx,1);
P(:,2) = reshape(Y,ny*nx,1);

% % use gptoolbox's in_element first output
% [I,~,~,~] = in_element([V(:,1) V(:,2)],F,P);
% % in_element is returning up to one triangle per point
% % (for some points there should be more than one triangle)

% tsurf(F,V,'FaceAlpha',0.1); axis equal;
% hold on

Z = zeros(ny*nx,1);
for k=1:ny*nx
    
    B = barycentric_coordinates(ones(size(F,1),1)*P(k,:),V(F(:,1),1:2),...
        V(F(:,2),1:2),V(F(:,3),1:2));
    
%     sum(I(k,:))
%     tris = find(I(k,:))

    tris = find((B(:,1)>=0)&(B(:,1)<=1)&(B(:,2)>=0)&(B(:,2)<=1)...
        &(B(:,3)>=0)&(B(:,3)<=1));
    
    if numel(tris)==0
        Z(k) = min_z;
    else
        max_height = -inf;
        for i=1:numel(tris)
            f = tris(i);
            curr_height = (V(F(f,1),3)+V(F(f,2),3)+V(F(f,3),3))/3;
            if curr_height>max_height
                max_height = curr_height;
            end
        end
        Z(k) = max_height;
    end
    
%     plot3(P(k,1),P(k,2),Z(k),'b.','MarkerSize',20)
    
end

sampled_HF = zeros(nx,ny,3);
sampled_HF(:,:,1) = X;
sampled_HF(:,:,2) = Y;
sampled_HF(:,:,3) = reshape(Z,nx,ny);
