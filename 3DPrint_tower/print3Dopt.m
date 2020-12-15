function [X,Vnew,F,fval,exitflag,output] = print3Dopt(V,F,x0,varargin)
%    This algorithm searches for the best orientation in
% space that solves the minimization problem:
%
%                   min sum W*||N-PrN||, 
%                   
% where N is the normal field of surface, PrN is the projection of the 
% normal field on the xy-plan and W is a weight function for each
% triangle of the mesh. This function uses fmincon to solve the minimiza-
% tion problem.
%
%Syntax:
%
% [X,Vnew,F,fval,exitflag,output] = print3Dopt(V,F,x0,varargin)
%
%Input:
%
%   V           #V by 3 vertex positions of the input mesh
%   F           #F by 3 vertex indices of each triangle of the input mesh
%   x0          initial point for fmincon
% Optional parameters:
%   zmin        minimum (fixed) height: it determines the height of the
%               supports
%   theta1_max  maximum theta1 angle for optimization (default is pi)
%   theta2_max  maximum theta2 angle for optimization (default is pi/2)
%
%Outputs:
%
%   X           2 by 1 vector solution of minimization problem
%   Vnew        #V by 3 list of vertices at the optimal position
%   F           #F by 3 list of triangle indices
%   fval        objective function value at the optimal point
%   exitflag    a value that describes the exit condition of fmincon
%   output      a structure output with information about the optimization 
%               process.

% It is better to read the file from outside the function
% and then call the function 
% [V,F]=readOBJ(surface); %% Reading .obj file from path

zmin = [];
theta1_max = pi;
theta2_max = pi/2;

v = 1;
while v<=numel(varargin) && ischar(varargin{v})
    switch varargin{v}
        case 'zmin'
            v = v+1;
            assert(v<=numel(varargin));
            zmin = varargin{v};
        case 'theta1_max'
            v = v+1;
            assert(v<=numel(varargin));
            theta1_max = varargin{v};
        case 'theta2_max'
            v = v+1;
            assert(v<=numel(varargin));
            theta2_max = varargin{v};
        otherwise
            break;
    end
    v = v+1;
end

% % % %Necessary parameters for fmincon
LB=[-theta1_max;-theta2_max]; %%Left box condition  %
UB=-LB; %% Right box condition        %
A=[];                                 %
B=[];                                 %
Aeq=[];                               %
Beq=[];                               %
% % % % % % % % % % % % % % % % % % % % 


Vh = rotatexy(V,F,[0;0],'center'); %% Translating the surface barycenter to
                                    % the origin.
N = normalsurf(Vh,F); %% Generating the normal field of the surface
Area = areatsurf(Vh,F); %% Generating a vector with triangle mesh areas

Bary = barycenter(Vh,F); 

% % % % % % % % % % % % % % % % % % %%fmincon applied to objective function
options = optimoptions('fmincon','Display','none');                       %
if isempty(zmin)
    [X,fval,exitflag,output] = fmincon(@(x)(...                               %
        (-sin(x(2,1)).*N(:,1)...                                              %
        + cos(x(2,1)).*sin(x(1,1)).*N(:,2)+...                                %
        cos(x(1,1)).*cos(x(2,1)).*N(:,3)>0).*...                              %
        Area.*(-sin(x(2,1))*...                                               %
        Bary(:,1)+ cos(x(2,1))*sin(x(1,1))*Bary(:,2) + ...                    %
          cos(x(1,1))*cos(x(2,1))*Bary(:,3) - min(-sin(x(2,1))...             %
          *Bary(:,1)+ cos(x(2,1))*sin(x(1,1))*Bary(:,2) + ...                 %
          cos(x(1,1))*cos(x(2,1))*Bary(:,3))))'*((sin(x(2,1)).^2)...          %
          .*(N(:,1).^2) + (cos(x(1,1)).^2)*(cos(x(2,1)).^2)...                %
          .*(N(:,3).^2) - 2.*sin(x(2,1)).*cos(x(2,1)).*sin(x(1,1))...         %
          .*N(:,1).*N(:,2) + (cos(x(2,1)).^2).*(sin(x(1,1)).^2)...            %
          .*(N(:,2).^2) + 2.*(cos(x(2,1)).^2).*sin(x(1,1))...                 %
          .*cos(x(1,1)).*N(:,2).*N(:,3) - 2.*sin(x(2,1))...                   %
          .*cos(x(2,1)).*cos(x(1,1)).*N(:,1)...                               %
          .*N(:,3)),x0,A,B,Aeq,Beq,LB,UB,[],options);                         %
else
    [X,fval,exitflag,output] = fmincon(@(x)(...                               %
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
          .*N(:,3)),x0,A,B,Aeq,Beq,LB,UB,[],options);                         %
end
    
  % This previous version uses the x,y,z coordinates of the first vertex
  % of each face, instead of x,y,z coodrinates of the barycenter
% options = optimoptions('fmincon','Display','none');                       %
% [X,fval,exitflag,output] = fmincon(@(x)(Area.*(-sin(x(2,1))*...           %
%     Vh(F(:,1),1)+ cos(x(2,1))*sin(x(1,1))*Vh(F(:,1),2) + ...              %
%       cos(x(1,1))*cos(x(2,1))*Vh(F(:,1),3) - min(-sin(x(2,1))...          %
%       *Vh(F(:,1),1)+ cos(x(2,1))*sin(x(1,1))*Vh(F(:,1),2) + ...           %
%       cos(x(1,1))*cos(x(2,1))*Vh(F(:,1),3))))'*((sin(x(2,1)).^2)...       %
%       .*(N(:,1).^2) + (cos(x(1,1)).^2)*(cos(x(2,1)).^2)...                %
%       .*(N(:,3).^2) - 2.*sin(x(2,1)).*cos(x(2,1)).*sin(x(1,1))...         %
%       .*N(:,1).*N(:,2) + (cos(x(2,1)).^2).*(sin(x(1,1)).^2)...            %
%       .*(N(:,2).^2) + 2.*(cos(x(2,1)).^2).*sin(x(1,1))...                 %
%       .*cos(x(1,1)).*N(:,2).*N(:,3) - 2.*sin(x(2,1))...                   %
%       .*cos(x(2,1)).*cos(x(1,1)).*N(:,1)...                               %
%       .*N(:,3)),x0,A,B,Aeq,Beq,LB,UB,[],options);                         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

Vnew = rotatexy(V,F,X,'center_back'); 

end