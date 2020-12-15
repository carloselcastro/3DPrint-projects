function Vnew = rotatexy(V,F,theta,option)
%Rotates a surface V around x and y axis by angle thetax,thetay
%
%Sintax:
%
%   Vnew = rotatexy(V,F,theta)
%
% Inputs:
%
%   V           #V by 3 matrix of surface's vertex coordinates.
%   F           #F by 3  matrix of indices of surface's triangle 
%               corners.
%   theta       2 by 1 vector with x and y rotation angle.
%   option      final position option of surface. 
%               Set 'plate' if the surface need to be positioned on printer
%               plate. Set 'center_back' to have the mesh translated such
%               that the barycenter of the result is equal to the original
%
% Output:
%
%   Vnew        #V x 3 matrix of surface's vertex coordinates after rotation.
%

% Making the input vector a column vector
if size(theta,2) > 1                    %
    theta = theta';                     %
end                                     %
% % % % % % % % % % % % % % % % % % % % %

thetax = theta(1,1); % Separating coordinates
thetay = theta(2,1); %

% % % % Translating the barycenter of the surface to the origin
Xbari = (V(F(:,1),1) + V(F(:,2),1) + V(F(:,3),1))/3;          %
Ybari = (V(F(:,1),2) + V(F(:,2),2) + V(F(:,3),2))/3;          %
Zbari = (V(F(:,1),3) + V(F(:,2),3) + V(F(:,3),3))/3;          %
X = sum(Xbari)/length(Xbari);                                 %
Y = sum(Ybari)/length(Ybari);                                 %
Z = sum(Zbari)/length(Zbari);                                 %
                                                              %
T=[1,0,0,-X;0,1,0,-Y;0,0,1,-Z;0,0,0,1];                       %
TV = T*[V ones(size(V,1),1)]';                                %
TV = TV';                                                     %                        
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

 % % % % % % % % % % % % % % % Setting and applying the rotation matrix
ROT1 = [cos(thetay) 0 sin(thetay)  0;...                              %
        0 1 0 0;...                                                   %
        -sin(thetay) 0 cos(thetay)  0;...                             %
        0 0 0 1];                                                     %
                                                                      %
ROT2 = [1 0 0 0;...                                                   %
        0 cos(thetax) -sin(thetax)  0;...                             %
        0 sin(thetax)  cos(thetax)  0;...                             %
        0 0 0 1];                                                     %
                                                                      %
ROT = ROT1*ROT2;                                                      %
Vnew = ROT*TV';  % Applying the rotation matrix                       %
if strcmp(option,'plate') == 1 % Positioning the surface after rotation
   T2=[1,0,0,0;0,1,0,0;0,0,1,-min(Vnew(3,:));0,0,0,1];                %
   Vnew = T2*Vnew;   
elseif strcmp(option,'center_back') == 1
    Vnew(1,:) = Vnew(1,:)+X;
    Vnew(2,:) = Vnew(2,:)+Y;
    Vnew(3,:) = Vnew(3,:)+Z;
end                                                                   %
Vnew = Vnew';      % New surface vertices matrix after rotation       %
Vnew = Vnew(:,1:3);%                                                  %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

end