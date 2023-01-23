function [dxb] = sys2(t,x,A,B,R,rT,SVPf,t1)
%function for calculating u(t) and x(t) and intergrating both using ode45

% initialize the solution output
dxb = zeros(4,1);  


s1 = interp1(t1,SVPf(:,1),t);
s2 = interp1(t1,SVPf(:,2),t);
s3 = interp1(t1,SVPf(:,3),t);
s4=interp1(t1,SVPf(:,4),t);
s5=interp1(t1,SVPf(:,5),t);
s6=interp1(t1,SVPf(:,6),t);
s7=interp1(t1,SVPf(:,7),t);
s8=interp1(t1,SVPf(:,8),t);
s9=interp1(t1,SVPf(:,9),t);
s10=interp1(t1,SVPf(:,10),t);
V1 = interp1(t1,SVPf(:,11),t);
V2 = interp1(t1,SVPf(:,12),t);
V3=interp1(t1,SVPf(:,13),t);
V4=interp1(t1,SVPf(:,14),t);
P = interp1(t1,SVPf(:,15),t);
                     
                     
% Define the S matrix
S = [s1 s2 s4 s5;
      s2 s3 s6 s7;
      s4 s6 s8 s9;
      s5 s7 s9 s10];

% Define the auxiliary vector v
V = [V1;
     V2;
     V3;
     V4];
 
% The states
x = [x(1);
     x(2);
     x(3);
     x(4)];
 
% The feedback gain
K = (R)^-1*B'*S;

% The control
u = -(K - (R)^-1*B'*V*(P)^-1*V')*x - (R)^-1*B'*V*(P)^-1*rT;

% The plant
dx = A*x + B*u;

dxb(1) = dx(1);
dxb(2) = dx(2);
dxb(3) = dx(3);
dxb(4) = dx(4);

end
