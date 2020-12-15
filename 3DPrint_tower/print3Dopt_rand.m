function [Xmin,Vnew,F,minim] = print3Dopt_rand(V,F,varargin)
%    This algorithm searches for the best orientation in
% space that solves the minimization problem:
%
%                   min sum W*||N-PrN||, 
%                   
% where N is the normal field of surface, PrN is the projection of the 
% normal field on the xy-plan and W is a weight function for each
% triangle of the mesh. This function uses fmincon to solve the minimiza-
% tion problem with a randomic multistart point strategy.
%
%Syntax:
%
% [Xmin,Vnew,F,minim] = print3Dopt_rand(V,F,varargin)
%
% Input:
%
%   V           #V by 3 vertex positions of the input mesh
%   F           #F by 3 vertex indices of each triangle of the input mesh
% Optional parameters:
%   zmin        minimum (fixed) height: it determines the height of the
%               supports
%   theta1_max  maximum theta1 angle for optimization (default is pi)
%   theta2_max  maximum theta2 angle for optimization (default is pi/2)
%
% Outputs:
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

m=32; %% Setting the number of random choices for initial point
% multistart = zeros(m,2); %%Creating a vector with random initial points

% % % % % % % % % % % % % % % % % % % % % % Begining of multistart strategy
for r=1:m
    
% % % % % % % % % % % % % % % % %Setting the random initial point
    ii=2*pi.*rand(1,1) -pi;                       %% x-coordinate
    jj=pi.*rand(1,1) - pi/2;                      %% y-coordinate                         
    x0=[ii;jj];  %% Setting the initial point with random choices

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
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 