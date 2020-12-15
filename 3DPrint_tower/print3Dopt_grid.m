function [Xmin,Vnew,F,minim] = print3Dopt_grid(V,F,varargin)
%    This algorithm searches for the best orientation in
% space that solves the minimization problem:
%
%                   min sum W*||N-PrN||, 
%                   
% where N is the normal field of surface, PrN is the projection of the 
% normal field on the xy-plan and W is a weight function for each
% triangle of the mesh. This function uses fmincon to solve the minimiza-
% tion problem with a uniform grid center multistart point strategy.
%
%Syntax:
%
% [Xmin,Vnew,F,minim] = print3Dopt_grid(V,F,varargin)
%
%Input:
%
%   V           #V by 3 vertex positions of the input mesh
%   F           #F by 3 vertex indices of each triangle of the input mesh
% Optional parameters:
%   zmin        minimum (fixed) height: it determines the height of the
%               supports
%   theta1_max  maximum theta1 angle for optimization (default is pi)
%   theta2_max  maximum theta2 angle for optimization (default is pi/2)
%
%Outputs:
%
%   Xmin        2 by 1 vector solution of the minimization problem
%   Vnew        #V by 3 list of vertices at the optimal position
%   F           #F by 3 list of triangle indices
%   minim       objective function value at the optimal point

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

minim=1e10; %% Setting a high minimum initial value
Xmin=zeros(2,1); %% Alocating the optimal point
h=pi/4; %% Setting the length to partition the grid domain

% % % % % % % % % % % % % % % % % % % % % % Begining of multistart strategy
for ii=-pi:h:pi-h                                                         
    for jj = -pi/2:h:pi/2-h                                               
        x0=[(ii+h)/2;(jj+h)/2]; %% Setting the initial point for fmincon  

        [X,Vnew_,F,fval,~,~] = print3Dopt(V,F,x0,'zmin',zmin,'theta1_max',...
        theta1_max, 'theta2_max', theta2_max);
        
% % %Getting the minimal value of function evaluation and the optimal point
        if fval < minim                                                   %
            minim = fval;                                                 %
            Xmin=X;                                                       %
            Vnew = Vnew_;                                                 %
        end                                                               %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  

    end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 