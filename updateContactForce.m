function [Ball1_new, Ball2_new] = updateContactForce(Ball1, Ball2)
global k

dist = sqrt((Ball1.pos_x - Ball2.pos_x)^2 + (Ball1.pos_y - Ball2.pos_y)^2);

if dist > Ball1.R + Ball2.R
    Ball1.F_c = 0;
    Ball1.n = [0;0];
    Ball2.F_c = 0;
    Ball2.n = [0;0];

else
    display('in contact')
    delta = Ball1.R + Ball2.R - dist;
    F_c = k*delta;
    N = [Ball1.pos_x-Ball2.pos_x; Ball1.pos_y-Ball2.pos_y];
    n = N/norm(N,2);
    Ball1.F_c = F_c;
    Ball1.n = n;
    Ball2.F_c = F_c;
    Ball2.n = -n;

end

Ball1_new = Ball1;
Ball2_new = Ball2;
