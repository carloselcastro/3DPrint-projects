function val = interpolate_height_field(sampled_HF,x,y)
%  This function interpolates a sampled height field at a point (x,y)
%
%Syntax:
%
% val = interpolate_height_field(sampled_HF,x,y)
%
%Input:
%
%   sampled_HF  nx by ny by 3 sampled height field where (:,:,1) stands for
%               x values, (:,:,2) stands for y values and (:,:,3) stands
%               for z (height) values
%   (x,y)       point location
%Outputs:
%  val          number that approximates sampled_HF(x,y)

min_x = min(min(sampled_HF(:,:,1)));
max_x = max(max(sampled_HF(:,:,1)));
min_y = min(min(sampled_HF(:,:,2)));
max_y = max(max(sampled_HF(:,:,2)));
min_z = min(min(sampled_HF(:,:,3)));

if (x>=min_x && x<=max_x && y>=min_y && y<=max_y)
    val = interp2(sampled_HF(:,:,1),sampled_HF(:,:,2),sampled_HF(:,:,3),x,y);
else
    val = min_z;
end