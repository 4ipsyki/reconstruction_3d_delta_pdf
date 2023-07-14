 function [c,ceq] = matrix_constraint(x)

c = [];

load('pars.mat','pars');
aa            = pars(18);
bb            = pars(19);
cc            = pars(20);
alpha         = pars(21);
beta          = pars(22);
gamma         = pars(23);

% axes 90 degrees apart
ceq7 = atan2d(norm(cross(x(:,2),x(:,3))),dot(x(:,2),x(:,3))) - alpha;
ceq8 = atan2d(norm(cross(x(:,3),x(:,1))),dot(x(:,3),x(:,1))) - beta;
ceq9 = atan2d(norm(cross(x(:,1),x(:,2))),dot(x(:,1),x(:,2))) - gamma;

ceq10 = norm(x(:,1))-norm(x(:,2)); % axes length equal
ceq11 = norm(x(:,2))-norm(x(:,3));
ceq12 = norm(x(:,3))-norm(x(:,1));

ceq13 = norm(x(:,1))-cc;
ceq14 = norm(x(:,2))-bb;

ceq = [ceq7, ceq8, ceq9, ceq10];

disp([norm(x(:,1)),norm(x(:,2)),norm(x(:,3))])
disp([ceq7+alpha,ceq8+beta,ceq9+gamma])

end
