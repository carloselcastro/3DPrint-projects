function R = print3Dopt_local_roations(V,F,theta_max)
% Given a mesh (result from the global rotation problem, 
% this code computes a rotation for each vertex that
% minimizes the overhang energy
%
% Example usage:
% R = print3Dopt_local_roations(V,F,theta_max)
%
% Input:
% V: #vertices by 3 list of vertex positions
% F: #faces by 3 list of face vertex indices
% theta_max: number  that determines the bounds for theta_1 
% and theta_2 (both belong to [-theta_max,theta_max])
%
% Output:
% R: #vertices by 3 b y 3 list of optimal per-vertex rotations
%


% % % %Necessary parameters for fmincon
LB=[-theta_max;-theta_max]; %%Left box condition  %
UB=-LB; %% Right box condition        %
A=[];                                 %
B=[];                                 %
Aeq=[];                               %
Beq=[];                               %
x0 = zeros(2,1);                      %
% % % % % % % % % % % % % % % % % % % % 

% fix zmin
zmin = min(V(:,3));
% fmincon options
options = optimoptions('fmincon','Display','none');  

% Initialize R as zeros
R = zeros(3,3,size(V,1));


% calculate optimal rotations for each vertex
for k=1:size(V,1)
    
    Flocal = mod(find(F==k)-1,size(F,1))+1;
    Flocal = F(Flocal,:);
    N = normalsurf(V,Flocal); %% Generating the normal field of the surface
    Area = areatsurf(V,Flocal); %% Generating a vector with triangle mesh areas
    Bary = barycenter(V,Flocal);
    
    [X,~,~,~] = fmincon(@(x)(...                               %
        (-sin(x(2,1)).*N(:,1)...                                              %
        + cos(x(2,1)).*sin(x(1,1)).*N(:,2)+...                                %
        cos(x(1,1)).*cos(x(2,1)).*N(:,3)>0).*...                              %
        Area.*(-sin(x(2,1))*...                                               %
        Bary(:,1)+ cos(x(2,1))*sin(x(1,1))*Bary(:,2) + ...                    %
        cos(x(1,1))*cos(x(2,1))*Bary(:,3) - zmin))'*((sin(x(2,1)).^2)...          %
        .*(N(:,1).^2) + (cos(x(1,1)).^2)*(cos(x(2,1)).^2)...                %
        .*(N(:,3).^2) - 2.*sin(x(2,1)).*cos(x(2,1)).*sin(x(1,1))...         %
        .*N(:,1).*N(:,2) + (cos(x(2,1)).^2).*(sin(x(1,1)).^2)...            %
        .*(N(:,2).^2) + 2.*(cos(x(2,1)).^2).*sin(x(1,1))...                 %
        .*cos(x(1,1)).*N(:,2).*N(:,3) - 2.*sin(x(2,1))...                   %
        .*cos(x(2,1)).*cos(x(1,1)).*N(:,1)...                               %
        .*N(:,3)),x0,A,B,Aeq,Beq,LB,UB,[],options);
    
    thetax = X(1);
    thetay = X(2);
    
    ROT1 = [cos(thetay) 0 sin(thetay); 0 1 0; -sin(thetay) 0 cos(thetay)];  
    ROT2 = [1 0 0; 0 cos(thetax) -sin(thetax); 0 sin(thetax) cos(thetax)]; 
    
    R(:,:,k) = ROT1*ROT2;
end