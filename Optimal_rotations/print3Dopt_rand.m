function [Xmin,Vnew,F,minim] = print3Dopt_rand(surface)
%    This algorithm uses the GPTOOLBOX (JACOBSON, 2018) readOBJ function 
% to read solids in .obj format and searches for the best orientation in
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
% [Xmin,V,F,minim] = print3Dopt_rand(surface)
%
%Input:
%
%   surface     path to .obj file
%
%Outputs:
%
%   Xmin        2 by 1 vector solution of the minimization problem
%   V           #V by 3 list of vertices at the optimal position
%   F           #F by 3 list of triangle indices
%   minim       objective function value at the optimal point
%
%JACOBSON, A. gptoolbox: Geometry processing toolbox. 
%http://github.com/alecjacobson/gptoolbox, 2018.
%

[V,F]=readOBJ(surface); %% Reading .obj file from path
minim=1e10; %% Setting a high minimum initial value
Xmin=zeros(2,1); %% Alocating the optimal point

m=32; %% Setting the number of random choices for initial point
% multistart = zeros(m,2); %%Creating a vector with random initial points

% % % %Necessary parameters for fmincon
LB=[-pi;-pi/2]; %%Left box condition  %
UB=-LB; %% Right box condition        %
A=[];                                 %
B=[];                                 %
Aeq=[];                               %
Beq=[];                               %
% % % % % % % % % % % % % % % % % % % % 

% % % % % % % % % % % % % % % % % % % % % % Begining of multistart strategy
for r=1:m
    
% % % % % % % % % % % % % % % % %Setting the random initial point
    ii=2*pi.*rand(1,1) -pi;                       %% x-coordinate
    jj=pi.*rand(1,1) - pi/2;                      %% y-coordinate                         
    x0=[ii;jj];%% Setting the initial point with random choices
%     multistart(r,:) = x0';%% Adding the initial point on vector
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    Vh = rotatexy(V,F,[0;0],'center');%% Translating the surface 
                                       % barycenter to the origin.
    N = normalsurf(Vh,F); %% Generating the normal field of the surface
    Area = areatsurf(Vh,F); %% Generating a vector with triangle mesh 
                             % areas
  
% % % % % % % % % % % % % % % % % % %%fmincon applied on objective function
    options = optimoptions('fmincon','Display','none');                   %
    [X,fval] = fmincon(@(x)(Area.*(-sin(x(2,1))*Vh(F(:,1),1)+ ...         %
        cos(x(2,1))*sin(x(1,1))*Vh(F(:,1),2) + ...                        %
        cos(x(1,1))*cos(x(2,1))*Vh(F(:,1),3) - ...                        %
        min(-sin(x(2,1))*Vh(F(:,1),1)+ ...                                %
        cos(x(2,1))*sin(x(1,1))*Vh(F(:,1),2) + ...                        %
        cos(x(1,1))*cos(x(2,1))*...                                       %
        Vh(F(:,1),3))))'*((sin(x(2,1)).^2).*(N(:,1).^2) + ...             %
        (cos(x(1,1)).^2)*(cos(x(2,1)).^2).*(N(:,3).^2) - ...              %
        2.*sin(x(2,1)).*cos(x(2,1)).*sin(x(1,1)).*N(:,1).*N(:,2) + ...    %
        (cos(x(2,1)).^2).*(sin(x(1,1)).^2).*(N(:,2).^2) + ...             %
        2.*(cos(x(2,1)).^2).*sin(x(1,1)).*cos(x(1,1)).*N(:,2).*N(:,3) - ...
        2.*sin(x(2,1)).*cos(x(2,1))...                                    %
        .*cos(x(1,1)).*N(:,1).*N(:,3)),x0,A,B,Aeq,Beq,LB,UB,[],options);  %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
 
% % %Getting the minimal value of function evaluation and the optimal point
        if fval < minim                                                   %
            minim = fval;                                                 %
            Xmin=X;                                                       %
        end                                                               %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
  
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

Vnew = rotatexy(V,F,Xmin,'plate');%% Rotating the surface by optimal angles
                                   % given by fmincon results, centering the
                                   % surface barycenter at the origin and
                                   % adjusting the surface to printer plate

function Vnew = rotatexy(V,F,theta,option)
%Rotate a surface V around x and y axis by angle thetax,thetay
%
% Vnew = rotatexy(V,F,theta)
% Inputs:
%   V #V x dim matrix of surface's vertex coordinates.
%   F #F x simplex_size  matrix of indices of surface's triangle corners.
%   theta 2x1 vector with x and y rotation angle.
% Output:
%   Vnew a #V x dim matrix of surface's vertex coordinates after rotation.
%Master thesis: Carlos 2019/2020
if size(theta,2) > 1
    theta = theta';
end

thetax = theta(1,1);
thetay = theta(2,1);

%putting surface in its barycenter
Xbari = (V(F(:,1),1) + V(F(:,2),1) + V(F(:,3),1))/3;
Ybari = (V(F(:,1),2) + V(F(:,2),2) + V(F(:,3),2))/3;
Zbari = (V(F(:,1),3) + V(F(:,2),3) + V(F(:,3),3))/3;
Xx = sum(Xbari)/length(Xbari);
Y = sum(Ybari)/length(Ybari);
Z = sum(Zbari)/length(Zbari);
% % % % % % % % % % % % % % % % % % % % % % % % % % 
T=[1,0,0,-Xx;0,1,0,-Y;0,0,1,-Z;0,0,0,1]; % Translation matrix if rotation
TV = T*[V ones(size(V,1),1)]';% % % % % % center to the origin.
TV = TV';

%rotate
ROT1 = [cos(thetay) 0 sin(thetay)  0;0 1 0 0;-sin(thetay) 0 cos(thetay)  0;0 0 0 1];
ROT2 = [1 0 0 0; 0 cos(thetax) -sin(thetax)  0;0 sin(thetax)  cos(thetax)  0;0 0 0 1];
ROT = ROT1*ROT2;             
Vnew = ROT*TV';  
if strcmp(option,'plate') == 1
   T2=[1,0,0,0;0,1,0,0;0,0,1,-min(Vnew(3,:));0,0,0,1];
   Vnew = T2*Vnew;
end

Vnew = Vnew';      % New V point after rotation
Vnew = Vnew(:,1:3);%
end %% Auxiliar function that 
                                               % rotates the surface.

function N = normalsurf(V,F)
%Find the Normal field of surface V.
%
% N = normalsurf(V,F)
% Inputs:
%   V #V x dim matrix of surface's vertex coordinates.
%   F #F x simplex_size  matrix of indices of surface's triangle corners.
% Output:
%   N #F x dim matrix of Normal field of surface.
%Master thesis: Carlos 2019/2020
v = V(F(:,3),:) - V(F(:,1),:);
w = V(F(:,2),:) - V(F(:,1),:);
vxw=cross(v,w);
N=normr(vxw);
end %% Auxiliar function that generate the
                                 % normal field of surface.

function A = areatsurf(V,F)
%Creates a vector with all triangle mesh areas of V.
%
% A = areatsurf(V,F)
% Inputs:
%   V #V x dim matrix of surface's vertex coordinates.
%   F #F x simplex_size  matrix of indices of surface's triangle corners.
% Output:
%   A #F x 1 vector of all triangle mesh areas of V.
%Master thesis: Carlos 2019/2020
a = sqrt((V(F(:,2),1) - V(F(:,1),1)).^2 + ...
         (V(F(:,2),2) - V(F(:,1),2)).^2 + ...
         (V(F(:,2),3) - V(F(:,1),3)).^2);
b = sqrt((V(F(:,3),1) - V(F(:,1),1)).^2 + ...
         (V(F(:,3),2) - V(F(:,1),2)).^2+ ...
         (V(F(:,3),3) - V(F(:,1),3)).^2);
c = sqrt((V(F(:,3),1) - V(F(:,2),1)).^2 + ...
         (V(F(:,3),2) - V(F(:,2),2)).^2+ ...
         (V(F(:,3),3) - V(F(:,2),3)).^2);
p = (a+b+c)./2;
            
A=sqrt(p.*(p-a).*(p-b).*(p-c));
end  %% Auxiliar function that generate a
                                 % vector with triangle mesh areas.

function [V,F,UV,TF,N,NF] = readOBJ(filename,varargin)
  % READOBJ reads an OBJ file with vertex/face information
  %
  % [V,F,UV,TF,N,NF] = readOBJ(filename)
  % [V,F,UV,TF,N,NF] = readOBJ(filename,'ParameterName',ParameterValue,...)
  %
  % Input:
  %  filename  path to .obj file
  %  Optional:
  %    'Quads' whether to output face information in X by 4 matrices (faces
  %      with degree larger than 4 are still triangulated). A trailing zero
  %      will mean a triangle was read.
  % Outputs:
  %  V  #V by 3 list of vertices
  %  F  #F by 3 list of triangle indices
  %  UV  #UV by 2 list of texture coordinates
  %  TF  #F by 3 list of triangle texture coordinates
  %  N  #N by 3 list of normals
  %  NF  #F by 3 list of triangle corner normal indices into N
  %
  % Examples:
  %   % read a quad/triangle mesh and display it
  %   [V,F] = readOBJ('quads.obj','Quads',true);
  %   % Turn triangles into degenerate quads 
  %   DF = (F==0).*F(:,[4 1 2 3])+(F~=0).*F;
  %   trisurf(DF,V(:,1),V(:,2),V(:,3));
  %
  %
  % See also: load_mesh, readOBJfast, readOFF

  % default values
  quads = false;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Quads'}, ...
    {'quads'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this workspace 
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end

  numv = 0;
  numf = 0;
  numuv = 0;
  numtf = 0;
  numn = 0;
  numnf = 0;

  % simplex size
  if quads
    ss = 4;
  else
    ss = 3;
  end
  % Amortized array allocation
  V = zeros(10000,3);
  F = zeros(10000,ss);
  UV = zeros(10000,3);
  TF = zeros(10000,ss);
  N = zeros(10000,3);
  NF = zeros(10000,ss);

  triangulated = false;
  all_ss = true;
  fp = fopen( filename, 'r' );
  type = fscanf( fp, '%s', 1 );
  count = 0;
  while strcmp( type, '' ) == 0
      line = fgets(fp);
      if strcmp( type, 'v' ) == 1
          v = sscanf( line, '%lf %lf %lf' );
          numv = numv+1;
          if(numv>size(V,1))
            V = cat(1,V,zeros(10000,3));
          end
          V(numv,:) = [v(1:3)'];
      elseif strcmp( type, 'vt')
          v = sscanf( line, '%f %f %f' );
          numuv = numuv+1;
          if size(UV,2)>2 && length(v) == 2
              UV = UV(:,1:2);
          end
          if(numuv>size(UV,1))
            UV = cat(1,UV,zeros(10000,length(v)));
          end
          UV(numuv,:) = [v'];
      elseif strcmp( type, 'vn')
          n = sscanf( line, '%f %f %f' );
          numn = numn+1;
          if(numn>size(N,1))
            N = cat(1,N,zeros(10000,3));
          end
          N(numn,:) = [n'];
      elseif strcmp( type, 'f' ) == 1
          [t, count] = sscanf(line,'%d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d');
          if (count>2)
              tf = t(2:3:end);
              nf = t(3:3:end);
              t = t(1:3:end);
          else
              [t, count] = sscanf(line, '%d/%d %d/%d %d/%d %d/%d %d/%d');
              if (count>2)
                  tf = t(2:2:end);
                  t = t(1:2:end);
                  nf = -ones(numel(t),1);
              else
                [t, count] = sscanf(line, '%d//%d %d//%d %d//%d %d//%d %d//%d');
                if (count>2)
                    nf = t(2:2:end);
                    t = t(1:2:end);
                    tf = -ones(numel(t),1);
                else
                    [t, count] = sscanf( line, '%d %d %d %d %d %d %d %d %d %d %d\n' );
                    if (count>2)
                      tf = -ones(numel(t),1);
                      nf = -ones(numel(t),1);
                    else
                      [t, count] = sscanf( line, '%d// %d// %d// %d// %d// %d// %d// %d// %d// %d// %d//\n' );
                      tf = -ones(numel(t),1);
                      nf = -ones(numel(t),1);
                    end
                end
              end
          end
          t = t + (t<0).*   (numv+1);
          tf = tf + (tf<0).*(numuv+1);
          nf = nf + (nf<0).*(numn+1);

          assert(numel(t) == numel(tf));
          assert(numel(t) == numel(nf));
          if numel(t) > ss
            if ~triangulated
              warning('Trivially triangulating high degree facets');
            end
            triangulated = true;
          end
          j = 2;
          i = 1;
          %Vt = V(t,:);
          %[~,A] = affine_fit(Vt);
          %VtA = Vt*A;
          %VtA0 = Vt*A;
          %[~,alpha] = curvature(VtA);
          %flip = -sign(sum(alpha));
          %E = [1:size(VtA,1);2:size(VtA,1) 1]';
          %[dV,dF] = triangle(VtA,E,[]);
          %if size(dF,1)>2
          %  tsurf(dF,dV);
          %  hold on;
          %  plot(VtA0([1:end 1],1),VtA0([1:end 1],2),'LineWidth',3);
          %  hold off
          %  pause
          %end
          while true
            if numel(t) > ss
              corners = [1 2 3];

              %plot(VtA0([1:end 1],1),VtA0([1:end 1],2));
              %hold on;
              %plot(VtA([1:3],1),VtA([1:3],2),'LineWidth',3);
              %hold off;
              %expand_axis(2);
              %pause;

              %[~,alpha] = curvature(VtA,[1 2;2 3]);
              %alpha = flip * alpha(2);
              %others = VtA(setdiff(1:end,corners),:);
              %these = VtA(corners,:);
              %w = inpolygon(others(:,1),others(:,2),these(:,1),these(:,2));
              %alpha
              %if alpha>=0 && ~any(w)
              %  % lazy
              %  t = t([2:end 1]);
              %  VtA = VtA([2:end 1],:);
              %  continue;
              %end
            else
              if all_ss && numel(t)<ss
                warning('Small degree facet found');
                all_ss = false;
              end
              corners = 1:numel(t);
            end
            numf = numf+1;
            if(numf>size(F,1))
              F = cat(1,F,zeros(10000,ss));
            end
            F(numf,1:numel(corners)) = [t(corners)'];
            numtf = numtf+1;
            if(numtf>size(TF,1))
              TF = cat(1,TF,zeros(10000,ss));
            end
            TF(numtf,1:numel(corners)) = [tf(corners)'];
            numnf = numnf+1;
            if(numnf>size(NF,1))
              NF = cat(1,NF,zeros(10000,ss));
            end
            NF(numnf,1:numel(corners)) = [nf(corners)'];
            if numel(t) <= ss
              break;
            end
            t = t([1 3:end]);
            tf = tf([1 3:end]);
            nf = nf([1 3:end]);
            %VtA = VtA([1 3:end],:);
            if numel(t) < 3
              break;
            end
          end
      elseif strcmp( type, '#' ) == 1
          %fscanf( fp, '%s\n', 1 );
          % ignore line
      end
      type = fscanf( fp, '%s', 1 );
  end
  fclose( fp );

  %try
  %    F = cell2mat(F);
  %end
  V = V(1:numv,:);
  F = F(1:numf,:);
  UV = UV(1:numuv,:);
  TF = TF(1:numtf,:);
  N = N(1:numn,:);
  NF = NF(1:numnf,:);

  %% transform into array if all faces have the same number of vertices

  if (size(UV,1)>0) UV = UV; end

end %% GPTOOLBOX
                                                           % function that
                                                           % read .obj
                                                           % files.
end