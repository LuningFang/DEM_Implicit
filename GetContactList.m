function Contacts = GetContactList(Balls, NB)
global k wall_L wall_R wall_B wall_T
nc = 1;
Contacts = [];

for i = 1:NB
    % make sure index i is always smaller than index j
    % contact detection with neightboring spheres
    for j = i+1:NB
        % apply contact detection
        dist = sqrt((Balls(i).pos_x - Balls(j).pos_x)^2 + (Balls(i).pos_y - Balls(j).pos_y)^2);
        delta = Balls(i).R + Balls(j).R - dist;
        if delta > 0
            Contacts(nc).index = nc;
            Contacts(nc).i = i;
            Contacts(nc).j = j;
            Contacts(nc).delta = delta;
            Contacts(nc).mag = k*delta;
            
            % contact force direction
            N = [Balls(i).pos_x-Balls(j).pos_x; Balls(i).pos_y-Balls(j).pos_y];
            n = N/norm(N,2);
            Contacts(nc).n1 = n(1);
            Contacts(nc).n2 = n(2);

            nc = nc + 1;
        end
    end
    
    % contact detection with boundaries
    if Balls(i).pos_x - wall_L < Balls(i).R
        Contacts(nc).index = nc;
        Contacts(nc).i = i;
        Contacts(nc).j = 0;
        Contacts(nc).n1 = 1;
        Contacts(nc).n2 = 0;
         Contacts(nc).delta = Balls(i).R - Balls(i).pos_x + wall_L;
         Contacts(nc).mag = k * Contacts(nc).delta;
        nc = nc + 1;
    end
    
    if -Balls(i).pos_x + wall_R < Balls(i).R
        Contacts(nc).index = nc;
        Contacts(nc).i = i;
        Contacts(nc).j = 0;
        Contacts(nc).n1 = -1;
        Contacts(nc).n2 = 0;
         Contacts(nc).delta = Balls(i).R + Balls(i).pos_x - wall_R;
         Contacts(nc).mag = k * Contacts(nc).delta;
        nc = nc + 1;
    end

    if Balls(i).pos_y - wall_B < Balls(i).R
        Contacts(nc).index = nc;
        Contacts(nc).i = i;
        Contacts(nc).j = 0;
        Contacts(nc).n1 = 0;
        Contacts(nc).n2 = 1;
         Contacts(nc).delta = (Balls(i).R - Balls(i).pos_y + wall_B);
         Contacts(nc).mag = k * Contacts(nc).delta;
        nc = nc + 1;
    end

    if -Balls(i).pos_y + wall_T < Balls(i).R
        Contacts(nc).index = nc;
        Contacts(nc).i = i;
        Contacts(nc).j = 0;
        Contacts(nc).n1 = 0;
        Contacts(nc).n2 = -1;
         Contacts(nc).delta = Balls(i).R + Balls(i).pos_y - wall_T;
         Contacts(nc).mag = k * Contacts(nc).delta;
        nc = nc + 1;
    end

    
    
end
    

end