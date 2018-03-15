function Phi = EvaluateRHS(Balls, Contacts)
global dt NB k wall_L wall_R wall_B wall_T
Phi = zeros(2*NB,1);

% without contact
for i = 1:NB
    Phi(2*i-1) = Phi(2*i-1) - Balls(i).Fx + Balls(i).m * Balls(i).acc_x;
    Phi(2*i) = Phi(2*i) - Balls(i).Fy+ Balls(i).m * Balls(i).acc_y;
end

% with contact
for nc = 1:length(Contacts)
    i = Contacts(nc).i;
    j = Contacts(nc).j;
    
    n1 = Contacts(nc).n1;
    n2 = Contacts(nc).n2;
    R1 = Balls(i).R;
    a1 = Balls(i).acc_x;
    a2 = Balls(i).acc_y;
    x1 = Balls(i).pos_x;
    x2 = Balls(i).pos_y;
    v1 = Balls(i).vel_x;
    v2 = Balls(i).vel_y;

    x1_new = x1 + dt * v1 + dt^2 * a1;
    x2_new = x2 + dt * v2 + dt^2 * a2;
    
    % determine if it's a boundary condition contact
    if j > 0

        R2 = Balls(j).R;
        a3 = Balls(j).acc_x;
        a4 = Balls(j).acc_y;

        x3 = Balls(j).pos_x;
        x4 = Balls(j).pos_y;

        v3 = Balls(j).vel_x;
        v4 = Balls(j).vel_y;
        
        x3_new = x3 + dt * v3 + dt^2 * a3;
        x4_new = x4 + dt * v4 + dt^2 * a4;
        
        
        dist = sqrt((x1_new - x3_new)^2 + (x2_new - x4_new)^2);
        delta = R1 + R2 - dist;
        F = k*delta;
        
        Phi(1 + 2*(i-1)) = Phi(1 + 2*(i-1)) - F*(x1_new-x3_new)/dist;
        Phi(2 + 2*(i-1)) = Phi(2 + 2*(i-1)) - F*(x2_new-x4_new)/dist;
        Phi(1 + 2*(j-1)) = Phi(1 + 2*(j-1)) + F*(x1_new-x3_new)/dist;
        Phi(2 + 2*(j-1)) = Phi(2 + 2*(j-1)) + F*(x2_new-x4_new)/dist;
        
    else
        if n1 == 1 && n2 ==0
            delta = R1 - (x1_new - wall_L);
            Phi(1 + 2*(i-1)) = Phi(1 + 2*(i-1)) - k*delta;
        end
        
        if n1 == -1 && n2 == 0
            delta = R1 - (wall_R - x1_new);
            Phi(1 + 2*(i-1)) = Phi(1 + 2*(i-1)) + k*delta;
        end
        
        if n1 == 0 && n2 == 1
            delta = R1 - (x2_new - wall_B);
            Phi(2 + 2*(i-1)) = Phi(2 + 2*(i-1)) - k*delta;
        end
           
        if n1 == 0 && n2 == -1
            delta = R1 - (wall_T - x2_new);
            Phi(2 + 2*(i-1)) = Phi(2 + 2*(i-1)) + k*delta;
        end
        
    end
end