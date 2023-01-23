function dSvb  = RiccatiV2(t,SVb,A,B,C,Q,R,ref,tr)

% initialize the solution output
dSvb = zeros(14,1);    % Chosen 5 because we need s1,s2,s3,v1,v2.


% Define the S matrix
S = [SVb(1),SVb(2),SVb(4),SVb(5);
    SVb(2),SVb(3),SVb(6),SVb(7);
    SVb(4),SVb(6),SVb(8),SVb(9);
    SVb(5),SVb(7),SVb(9),SVb(10)];

% Define the auxiliary vector V
v = [SVb(11);
     SVb(12);
     SVb(13);
     SVb(14)];

% The feedback gain
K = (inv(R)*B'*S);


% Interpolate the reference trajectory
r = interp1(tr,ref,t); %tr = ref time, t = current time, r(tol)

% The Riccati equation
dS = (A'*S + S*A - S*B*inv(R)*B'*S + C'*Q*C);

% Auxiliary vector ODE
dv = (A - B*K)'*v + C'*Q*r;

dSvb(1)=dS(1,1);
dSvb(2)=dS(1,2);
dSvb(3)=dS(2,2);
dSvb(4)=dS(1,3);
dSvb(5)=dS(1,4);
dSvb(6)=dS(2,3);
dSvb(7)=dS(2,4);
dSvb(8)=dS(3,3);
dSvb(9)=dS(3,4);
dSvb(10)=dS(4,4);
dSvb(11)=dv(1,1);
dSvb(12)=dv(2,1);
dSvb(13)=dv(3,1);
dSvb(14)=dv(4,1);

end