% check Jacobian
global NB
epsilon = 1e-4;

phi_q_approx = zeros(4*NB, 4*NB);
for i =1:4*NB
   dq = zeros(4*NB,1);
   dq(i) = epsilon;
   phi = EvaluateRHS(Balls, Contacts, prev_val, new_val);
   phi_new = EvaluateRHS(Balls, Contacts, prev_val, new_val+dq);
   phi_q = (phi_new - phi)/epsilon;
   phi_q_approx(:,i) = phi_q;
end

phi_q_approx - Phi_q