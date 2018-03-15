function newBall = updatePosition(Ball)
global dt
newBall = Ball;
newBall.pos_x = Ball.pos_x + dt*Ball.vel_x;
newBall.pos_y = Ball.pos_y + dt*Ball.vel_y;
end