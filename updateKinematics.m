function    Balls_new = updateKinematics(Balls, new_val)
global NB
Balls_new = Balls;

for i = 1:NB
    Balls_new(i).acc_x = new_val(2*i - 1);
    Balls_new(i).acc_y = new_val(2*i - 0);
end