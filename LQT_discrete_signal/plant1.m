function dxs = plant1(t,x1,A,B,R,SVf,t1)

% initialize the solution output
dxs = zeros(4,1);  


s1 = interp1(t1,SVf(:,1),t);
s2 = interp1(t1,SVf(:,2),t);
s3 = interp1(t1,SVf(:,3),t);
s4=interp1(t1,SVf(:,4),t);
s5=interp1(t1,SVf(:,5),t);
s6=interp1(t1,SVf(:,6),t);
s7=interp1(t1,SVf(:,7),t);
s8=interp1(t1,SVf(:,8),t);
s9=interp1(t1,SVf(:,9),t);
s10=interp1(t1,SVf(:,10),t);
V1 = interp1(t1,SVf(:,11),t);
V2 = interp1(t1,SVf(:,12),t);
V3=interp1(t1,SVf(:,13),t);
V4=interp1(t1,SVf(:,14),t);


% Define the S matrix
S = [s1 s2 s4 s5;
      s2 s3 s6 s7;
      s4 s6 s8 s9;
      s5 s7 s9 s10];

% Define the auxiliary vector v
v = [V1;
     V2;
     V3;
     V4];
 
% The feedback gain
K = R^-1*B'*S;


% The states
x = [x1(1);
     x1(2);
     x1(3);
     x1(4)];

% The control
u = -K*x + R^-1*B'*v;

% The plant
dx = A*x + B*u;


dxs(1) = dx(1);
dxs(2) = dx(2);
dxs(3) = dx(3);
dxs(4) = dx(4);


end