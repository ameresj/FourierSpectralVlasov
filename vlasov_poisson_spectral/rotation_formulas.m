syms theta

R(theta)=[ cos(theta), -sin(theta); sin(theta) , cos(theta)]

pi=sym('pi');

% check antiderivative
simplify(diff(R(theta-pi/2),theta)-R(theta))

% forward scovel in x
simplify([R(theta - pi/2)- R(0 -pi/2)]*R(-theta) )


% adjoint scovel in x


R(0)

simplify(R(0 -pi/2))