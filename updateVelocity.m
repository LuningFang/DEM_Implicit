function newBall = updateVelocity(Ball)
global dt
vel = [Ball.vel_x; Ball.vel_y] + dt*Ball.F_c*Ball.n/Ball.m;
newBall = Ball;
newBall.vel_x = vel(1);
newBall.vel_y = vel(2);
end