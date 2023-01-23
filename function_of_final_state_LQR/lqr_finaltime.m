%% 

clc
clear
close all

%% plant

b=340.3;
B=[0,b,0,0]';
wn=3.3;
zeta=0.093;
k=4.748;
c=1.168;
A=[0,1,0,0;
    -wn^2,-2*zeta*wn,0,0;
    0,0,0,1;
    -k,0,0,-c];

C=[0 0 1 0];%output


%% cost index parameters

%Transition state weight
Q=(100)*[1,0,0,0;
    0,1,0,0;
    0,0,1,0;
    0,0,0,1]; %100

% Control weight
R=1e7;%

%Final state weight 
ST=1*[1,0,0,0;
    0,1,0,0;
    0,0,1,0;
    0,0,0,1];
%final time
T=15;

%initial condition
x0=[0,0,0,0]';
%final modified state transition matrix
VT = C';

% Reachability gramian at final time
PT = 1e-4;%  % Set very small to avoid singularity

%% Solve the system to get S(t), V(t), and P(t)

% Solve the ODEs backwards in time 
[t,SVP1] = ode45(@(t,SVP) sys1(t,SVP,A,B,R,Q),[0 T],[ST(1,1) ST(1,2) ST(2,2) ST(1,3) ST(1,4) ST(2,3) ST(2,4) ST(3,3) ST(3,4) ST(4,4) VT(1) VT(2) VT(3) VT(4) PT]);
%t is time vector and SVP contains information for S(t), V(t), and P(t)

% Flip the solution to have S(t),V(t),and P(t) forward in time
N = length(t);

%makes matrices filled with zeros 
Kt = zeros(N,4);
SVPf = zeros(N,15); 
t1 = zeros(1,N);

for k = 1:N
    
    S1= [SVP1(N-k+1,1),SVP1(N-k+1,2),SVP1(N-k+1,4),SVP1(N-k+1,5);
              SVP1(N-k+1,2),SVP1(N-k+1,3),SVP1(N-k+1,6),SVP1(N-k+1,7);
               SVP1(N-k+1,4),SVP1(N-k+1,6),SVP1(N-k+1,8),SVP1(N-k+1,9);
               SVP1(N-k+1,5),SVP1(N-k+1,7),SVP1(N-k+1,9),SVP1(N-k+1,10)]; 
           
    Kt(k,:) = inv(R)*B'*S1; % calculate K forward in time
    
    SVPf(k,:) = SVP1(N-k+1,:); %flips matrix SVP to get S(t), V(t), and P(t) forward in time. 
                              %SVPf contains information for S(t), V(t), and P(t)
    t1(k) = T - t(N-k+1); %flips time interval 
end


%% Plot S, V, P, K

figure(1)
plot(t1,SVPf(:,1:10),'linewidth',1.5)
grid on
xlabel('Time (sec)')
ylabel('S')
legend('s_1','s_2','s_3','s_4','s_5','s_6','s_7','s_8','s_9','s_{10}','location','northeastoutside')
% title('Kernal matrix trajectory S(t)')


figure(2)
plot(t1,SVPf(:,11:14),'linewidth',1.5)
grid on
xlabel('Time (sec)')
ylabel('V')
legend('V1(t)','V2(t)','V3(t)','V4(t)','location','northwest')
% title('Modified state transition matrix trajectory V(t)')

figure(3)
plot(t1,SVPf(:,15),'linewidth',1.5)
grid on
xlabel('Time (sec)')
ylabel('P')
% title('Reachability gramian P(t)')

figure(4)
plot(t1,Kt,'linewidth',1.5)
grid on
xlabel('Time (sec)')
ylabel('K')
legend('K1(t)','K2(t)','K3(t)','K4(t)','location','northwest')
% title('Feedback Gain K(t)')

%%  Using MATLAB simulation to find the optimal tracking control 
%to track the reference r(t)=sin(t), please plot the control u(t) and x(t), 
%r(t) versus time.

%refrence signal (r(t)=t^2) at final time T
r=0.25*sin(t);
rT=r(end);



%%  Since we have the values V(t), S(t), and P(t), we can solve the plant forward in time.


% Solve the tracking problem forward in time to get x(t) and u(t)

[tcl,xf] = ode45(@(t,x) sys2(t,x,A,B,R,rT,SVPf,t1),t,[x0(1) x0(2) x0(3) x0(4)]); 
%tcl is time and xuf contains information for x(t) and u(t)

%%  to integrate u(t) forward in time, we calculate u(t) initial
%from initial values of: S(t), V(t), K(t), P(t) and x(t)
for k=1:N-1
    
S = [SVPf(k,1),SVPf(k,2),SVPf(k,4),SVPf(k,5);
      SVPf(k,2),SVPf(k,3),SVPf(k,6),SVPf(k,7);
      SVPf(k,4),SVPf(k,6),SVPf(k,8),SVPf(k,9);
      SVPf(k,5),SVPf(k,7),SVPf(k,9),SVPf(k,10)];  

K = inv(R)*B'*S;

V=[SVPf(k,11);SVPf(k,12);SVPf(k,13);SVPf(k,14)]; 
P=SVPf(k,15);

%calculating initial value of u(t)

u0 = -(K - (R)^-1*B'*V*(P)^-1*V')*xf(k+1,:)' - (R)^-1*B'*V*(P)^-1*rT;
u(k)=u0;
end

%%  Plot closed loop trajectory x(t), r(t) and u(t)
figure(5)

plot(tcl,xf(:,1:4),'linewidth',2)
hold on 
plot(tcl,0.25*sin(tcl),'linewidth',2)
grid on
legend('\theta','q','x','vx','r(t)','fontsize',15,'location','southwest')
% title('Closedloop State Trajectory','fontsize',15)
xlabel('Time(s)','fontsize',15)
ylabel('X','fontsize',15)
figure(6)
plot(tcl(1:end-1),u,'linewidth',2)
grid on
% title('Optimal control','fontsize',15)
xlabel('Time(s)','fontsize',15)
ylabel('u(t)','fontsize',15)


