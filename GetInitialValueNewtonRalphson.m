function initial_val = GetInitialValueNewtonRalphson(Balls)

NB = length(Balls);
initial_val = zeros(2*NB,1);

for i = 1:NB
    initial_val(0*NB + 1 + 2*(i-1)) = Balls(i).acc_x;
    initial_val(0*NB + 2 + 2*(i-1)) = Balls(i).acc_y;
end