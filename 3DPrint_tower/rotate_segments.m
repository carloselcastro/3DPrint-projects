function segs = rotate_segments(V,F,theta1_max,theta2_max,seg)
%   Given a mesh and segmentation, finds the optimal rotation for each
%   segment within the angle range 
%   [-theta1_max,theta1_max] x [-theta2_max,theta2_max]
%
% Example usage:
%       rotated_segs = rotate_segments(V,F,pi/6,pi/12,seg);
%
% Input:
%
%   V           #V by 3 vertex positions of the input mesh
%   F           #F by 3 vertex indices of each triangle of the input mesh
%   theta1_max  maximum theta1 angle for optimization (default is pi)
%   theta2_max  maximum theta2 angle for optimization (default is pi/2)
%   seg         vector containing segment index for each face
% 
% Output:
%
%   segs: vertex coordinates of each segment

n_segs = max(seg);
colors = rand(n_segs,3);

for k=1:n_segs
    
    fprintf('rotating segment %d of %d \n', k,n_segs)
    
    [~,Vnew,~,~] = print3Dopt_grid(V,F(seg==k,:),'zmin',...
        min(V(:,3)),'theta1_max',theta1_max,'theta2_max',theta2_max);
    
    segs{k} = Vnew;
    tsurf(F(seg==k,:),Vnew,'FaceColor',colors(k,:)); axis equal; hold on
    view([1 0 0])
        
end