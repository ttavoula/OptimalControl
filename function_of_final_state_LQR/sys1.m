function dsol  = sys1(~,SVP,A,B,R,Q)
%ODEs for S, V, and P

% initialize the solution output
dsol=zeros(15,1);        % Chosen 15 because we need s1->s10,V1->V4, and P


% Define the S matrix
S = [SVP(1),SVP(2),SVP(4),SVP(5);
    SVP(2),SVP(3),SVP(6),SVP(7);
    SVP(4),SVP(6),SVP(8),SVP(9);
    SVP(5),SVP(7),SVP(9),SVP(10)];

% Define the auxiliary vector V
V = [SVP(11);
     SVP(12);
     SVP(13);
     SVP(14)];

% The Riccati equation
dS=(A'*S + S*A - S*B*inv(R)*B'*S + Q);


dsol(1)=dS(1,1);
dsol(2)=dS(1,2);
dsol(3)=dS(2,2);
dsol(4)=dS(1,3);
dsol(5)=dS(1,4);
dsol(6)=dS(2,3);
dsol(7)=dS(2,4);
dsol(8)=dS(3,3);
dsol(9)=dS(3,4);
dsol(10)=dS(4,4);

% Feedback gain
K = (inv(R)*B'*S);

% Auxiliary vector ODE
dV=((A - B*K)'*V);

dsol(11)=dV(1,1);
dsol(12)=dV(2,1);
dsol(13)=dV(3,1);
dsol(14)=dV(4,1);

% The gramian
dP= -1*(V'*B*inv(R)*B'*V);%for dP we put -ve sign

dsol(15)=dP;

end