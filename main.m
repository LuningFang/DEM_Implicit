clc
clear all

global k dt NB Radius wall_L wall_R wall_B wall_T
tic;
% system parameters
k = 1e9;
NB = 4;
Mass = 1*ones(NB,1);
Radius = 1*ones(NB,1);
wall_L = 0;
wall_R = 4;
wall_B = 0;
wall_T = 4;
%F_ext_array = [0;0;0;0;0;0;-20;-20];


% simulation parameters
Tend = 0.005;
dt = 1e-5;


for i = 1:NB
    Balls(i).index = i;
    Balls(i).pos_x = 0;
    Balls(i).pos_y = 0;
    Balls(i).vel_x = 0;
    Balls(i).vel_y = 0;
    Balls(i).acc_x = 0;
    Balls(i).acc_y = 0;
    Balls(i).R = Radius(i);
    Balls(i).m = Mass(i);
    Balls(i).Fx = 0;
    Balls(i).Fy = 0;
end

time = 0:dt:Tend;

Balls(1).pos_x = 1;
Balls(1).pos_y = 1;
Balls(2).pos_x = 1;
Balls(2).pos_y = 3;
Balls(3).pos_x = 3;
Balls(3).pos_y = 1;
Balls(4).pos_x = 3;
Balls(4).pos_y = 3;
Balls(4).vel_x = 0;
Balls(4).vel_y = 0;
Balls(4).Fx = -100;
Balls(4).Fy = -100;


RenderInfo = zeros(length(time), 2*NB);
total_delta = zeros(length(time),1);


stateVar = GetInitialValueNewtonRalphson(Balls);
RenderInfo(1,:) = stateVar(1:2*NB);
total_delta(1) = 0;

for i = 2:length(time)
    Contacts = GetContactList(Balls, NB);
    
    
    prev_val = GetInitialValueNewtonRalphson(Balls);
    %tol = 1e-13/(dt^2);
    tol = 1e-3;
    err = 10e4;
    new_val = prev_val; % [x1 ... xn y1 ... yn] at previous time step
    num_itr = 1;
    while err > tol
        Phi_q = AssembleJacobian(Balls, Contacts);
        
        Phi = EvaluateRHS(Balls, Contacts);
        new_val = new_val - inv(Phi_q)*Phi;
        Balls_new = updateKinematics(Balls, new_val);
        Balls = Balls_new;

        Phi_new = EvaluateRHS(Balls, Contacts);
        err = norm(Phi_new,2);
        
        % check if the matrix is diagonal dominant or not
%         isDiagonalDominant = 0;
%         if (IsDiagonalDominant(Phi_q) == 1)
%             isDiagonalDominant = 1;
%         end
        
%         isPositiveDefinite = 0;
%         [~,p] = chol(Phi_q);
%         if p == 0
%             isPositiveDefinite = 1;
%         end
%         fprintf('time = %f, err = %g, diagonalDominant = %d, posDef = %d, max_element = %f\n', time(i), err, isDiagonalDominant, isPositiveDefinite, max(max(Phi_q)));
        
        num_itr = num_itr + 1;
        if num_itr > 20
            break;
        end
    end
    
    Balls_new = updateKinematics(Balls, new_val);
    Balls = Balls_new;
    for nb = 1:NB
        Balls(nb).vel_x = Balls(nb).vel_x + dt * Balls(nb).acc_x;
        Balls(nb).vel_y = Balls(nb).vel_y + dt * Balls(nb).acc_y;
        Balls(nb).pos_x = Balls(nb).pos_x + dt * Balls(nb).vel_x;
        Balls(nb).pos_y = Balls(nb).pos_y + dt * Balls(nb).vel_y;
        RenderInfo(i, 2*nb-1) = Balls(nb).pos_x;
        RenderInfo(i, 2*nb  ) = Balls(nb).pos_y;
    end
    total_delta(i) = getTotalPenetration(Contacts);
    

    
end

computeTime = toc;
figure;
plot(time,total_delta);
xlabel('time(sec)')
ylabel(strcat('total penetration of all contacts, ', sprintf('\\delta', 'Interpreter', 'latex')))

title(sprintf('$$\\Delta t = %g, k = %g, \\delta_{ss} = %g$$', dt, k, total_delta(end)), 'Interpreter', 'latex');
str_dir = pwd;
str_figname = sprintf('/dt_1e%d_Implicit.png', log10(dt));
xlim = max(get(gca,'XLim'));
ylim = max(get(gca,'YLim'));
text(0.6*xlim, 0.95*ylim, 'Implicit Scheme');
text(0.6*xlim, 0.9*ylim, sprintf('TimeCost = %.2fsec', computeTime));

print(gcf, strcat(str_dir, str_figname), '-dpng', '-r300');


