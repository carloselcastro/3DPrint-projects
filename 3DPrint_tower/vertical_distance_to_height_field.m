function D = vertical_distance_to_height_field(P,sampled_HF)
%  This function approximates the vertical distance from points
%  to a given height field. If the point is inside the graph of the
%  the height field, then it returns infinity
%
% Syntax:
%
% D = vertical_distance_to_height_field(P,sampled_HF)
%
%Input:
%
%   P           #P by 3 (x,y,z) coordinates of query points
%   sampled_HF  nx by ny by 3 sampled height field where (:,:,1) stands for
%               x values, (:,:,2) stands for y values and (:,:,3) stands
%               for z (height) values
%Outputs:
%  D          #P by 1 vector with vertical distances to query points

% TO-DO: vectorize this code (will probably have to vectorize 
% interpolate_height_field too

D = zeros(size(P,1),1);
for k=1:size(P,1)
    % approximate height field value at (x,y)-coordinates of the point
    z_HF = interpolate_height_field(sampled_HF,P(k,1),P(k,2));
    % if height field is lower return the difference
    if (z_HF <= P(k,3))
        D(k) = P(k,3)-z_HF;
    % else return infinity (point P is inside the height field
    else
        D(k) = inf;
    end
end