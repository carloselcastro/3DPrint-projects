function V_all = print3Dopt_local_global(V,F,theta_max,dec_factor)

V_prev = zeros(size(V));
V_all(:,:,1) = V;
minim = inf;

while (size(V_all,3)==1 || minim<dec_factor*minim_prev)
    
    minim_prev = minim;
    V_prev = V;
    % global step
    [X,V_global,~,minim] = print3Dopt_grid(V,F);
    fprintf('optimal angles: %.4f %.4f \n',X(1),X(2))
    fprintf('energy after global step: %.6f \n',minim)
    V_all(:,:,end+1) = V_global;
    
    % local step
    R = print3Dopt_local_roations(V_global,F,theta_max);
    V = arap_positions(V_global,F,R);
    V_all(:,:,end+1) = V;
    
    fprintf('end of iteration %d \n',(size(V_all,3)-1)/2)
%     figure; tsurf(F,V); axis equal;
%     input('');
    
end