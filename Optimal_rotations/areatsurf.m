function A = areatsurf(V,F)
%Creates a vector with all triangle mesh areas of surface.
%
%Sintax:
%
%   A = areatsurf(V,F)
%
% Inputs:
%
%   V           #V x 3 matrix of surface's vertex coordinates.
%   F           #F x 3  matrix of indices of surface's triangle corners.
%
% Output:
%
%   A           #F x 1 vector of all triangle mesh areas of V.
%

% Calculating parameters for triangle area using Heron's formula
a = sqrt((V(F(:,2),1) - V(F(:,1),1)).^2 + ...  %triangle edge a 
         (V(F(:,2),2) - V(F(:,1),2)).^2 + ...                 %
         (V(F(:,2),3) - V(F(:,1),3)).^2);                     %
b = sqrt((V(F(:,3),1) - V(F(:,1),1)).^2 + ...  %triangle edge b
         (V(F(:,3),2) - V(F(:,1),2)).^2+ ...                  %
         (V(F(:,3),3) - V(F(:,1),3)).^2);                     %
c = sqrt((V(F(:,3),1) - V(F(:,2),1)).^2 + ...  %triangle edge c
         (V(F(:,3),2) - V(F(:,2),2)).^2+ ...                  %
         (V(F(:,3),3) - V(F(:,2),3)).^2);                     %
p = (a+b+c)./2; %Semiperimeter                                %
                                                              %
A=sqrt(p.*(p-a).*(p-b).*(p-c)); % Calculating area of triangles
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

end