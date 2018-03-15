function Phi_q = AssembleJacobian(Balls, Contacts)
global NB dt k

Phi_q = zeros(2*NB, 2*NB);

% assemble non-contact components
for i = 1:NB
    Phi_q(2*i-1,2*i-1) = Balls(i).m;
    Phi_q(2*i, 2*i) = Balls(i).m;
end


for nc = 1:length(Contacts)
    i = Contacts(nc).i;
    j = Contacts(nc).j;
    n1 = Contacts(nc).n1;
    n2 = Contacts(nc).n2;
    
    
    % set contact normal direction to be fixed, invariant wrt the variables
    %     dFdx = k/sqrt((x3-x1)^2 + (x4-x2)^2)*[x3-x1, x4-x2, x1-x3, x2-x4];
    %     nmdt = [-n1/m1, -n2/m1, n1/m2, n2/m2]*dt;
    %     DFDX = dFdx'*nmdt;
    
    % contact normal dependent on displacement
    
    if j > 0
        m1 = Balls(i).m;
        m2 = Balls(j).m;
        R1 = Balls(i).R;
        R2 = Balls(j).R;
        a1 = Balls(i).acc_x;
        a2 = Balls(i).acc_y;
        a3 = Balls(j).acc_x;
        a4 = Balls(j).acc_y;

        x1 = Balls(i).pos_x;
        x2 = Balls(i).pos_y;
        x3 = Balls(j).pos_x;
        x4 = Balls(j).pos_y;

        v1 = Balls(i).vel_x;
        v2 = Balls(i).vel_y;
        v3 = Balls(j).vel_x;
        v4 = Balls(j).vel_y;
        

        % delta = distance in between the center of mass of spheres
        x1_new = x1 + dt * v1 + dt^2 *a1;
        x2_new = x2 + dt * v2 + dt^2 *a2;
        x3_new = x3 + dt * v3 + dt^2 *a3;
        x4_new = x4 + dt * v4 + dt^2 *a4;
        dist = sqrt((x1_new-x3_new)^2 + (x2_new-x4_new)^2);
        
        % derivative of inv of delta w.r.t. [x1, x2, x3, x4]
        d_inv_dist_dx1 = -1/dist^3*(x1_new-x3_new);
        d_inv_dist_dx2 = -1/dist^3*(x2_new-x4_new);
        d_inv_dist_dx3 =  1/dist^3*(x1_new-x3_new);
        d_inv_dist_dx4 =  1/dist^3*(x2_new-x4_new);
        
        % derivative of contact force in x and y direction, Fx and Fy
        % with respect to [a1, a2, a3, a4]
        d_Fx_da1 = dt^2 * (k*(R1+R2)*d_inv_dist_dx1*(x1_new-x3_new) + k*((R1+R2)/dist - 1));
        d_Fx_da2 = dt^2 * (k*(R1+R2)*d_inv_dist_dx2*(x1_new-x3_new));
        d_Fx_da3 = dt^2 * (k*(R1+R2)*d_inv_dist_dx3*(x1_new-x3_new) - k*((R1+R2)/dist - 1));
        d_Fx_da4 = dt^2 * (k*(R1+R2)*d_inv_dist_dx4*(x1_new-x3_new));
        
        d_Fy_da1 = dt^2 * (k*(R1+R2)*d_inv_dist_dx1*(x2_new-x4_new));
        d_Fy_da2 = dt^2 * (k*(R1+R2)*d_inv_dist_dx2*(x2_new-x4_new) + k*((R1+R2)/dist - 1));
        d_Fy_da3 = dt^2 * (k*(R1+R2)*d_inv_dist_dx3*(x2_new-x4_new));
        d_Fy_da4 = dt^2 * (k*(R1+R2)*d_inv_dist_dx4*(x2_new-x4_new) - k*((R1+R2)/dist - 1));
        
        DFDX = zeros(4,4);
        DFDX(1,:) = -[d_Fx_da1 d_Fx_da2 d_Fx_da3 d_Fx_da4];
        DFDX(2,:) = -[d_Fy_da1 d_Fy_da2 d_Fy_da3 d_Fy_da4];
        DFDX(3,:) =  [d_Fx_da1 d_Fx_da2 d_Fx_da3 d_Fx_da4];
        DFDX(4,:) =  [d_Fy_da1 d_Fy_da2 d_Fy_da3 d_Fy_da4];
        
        Phi_q(1+(i-1)*2:2+(i-1)*2, 1+(i-1)*2:2+(i-1)*2) = Phi_q(1+(i-1)*2:2+(i-1)*2, 1+(i-1)*2:2+(i-1)*2) + DFDX(1:2,1:2);
        Phi_q(1+(i-1)*2:2+(i-1)*2, 1+(j-1)*2:2+(j-1)*2) = Phi_q(1+(i-1)*2:2+(i-1)*2, 1+(j-1)*2:2+(j-1)*2) + DFDX(1:2,3:4);
        Phi_q(1+(j-1)*2:2+(j-1)*2, 1+(i-1)*2:2+(i-1)*2) = Phi_q(1+(j-1)*2:2+(j-1)*2, 1+(i-1)*2:2+(i-1)*2) + DFDX(3:4,1:2);
        Phi_q(1+(j-1)*2:2+(j-1)*2, 1+(j-1)*2:2+(j-1)*2) = Phi_q(1+(j-1)*2:2+(j-1)*2, 1+(j-1)*2:2+(j-1)*2) + DFDX(3:4,3:4);
        
    else
        if j == 0
            m1 = Balls(i).m;
            
            if n1 ~= 0 && n2 ==0
                
                Phi_q(1+(i-1)*2, 1+(i-1)*2) = Phi_q(1+(i-1)*2, 1+(i-1)*2) + k*dt^2;
            else
                if n1 == 0 && n2 ~= 0
                    Phi_q(2+(i-1)*2, 2+(i-1)*2) = Phi_q(2+(i-1)*2, 2+(i-1)*2) + k*dt^2;

                    
                end
            end
            
        end
    end
    
end