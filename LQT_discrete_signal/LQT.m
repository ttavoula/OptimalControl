%% This code is for simulation of LQT with discrete reference signal.
clc
clear all
close all

%% The Plant
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

% The output 
C=[0 0 1 0];
%% Performance index
% Transition state weight
Q = 10  ;   
% Control weight
R = 200; 
% Final state weight
P = 1; 
% Initial state
x0 = zeros(4,1);

%%
T = 15;     % Final time 

%% Reference trajectory

x0_val = [0, -0.25, 0.25, -0.25, 0, 0];

dt = .1;
T_period = 3;  % How long it is constant

N = length(x0_val);

t_vec = (0:N-1) * T_period;

figure
stairs(t_vec, x0_val)
grid

t_vec_over = t_vec(1):dt:t_vec(end);

x0_val_over = interp1(t_vec, x0_val, t_vec_over, 'previous');
hold on
stairs(t_vec_over, x0_val_over)

tr = t_vec_over;
r1 = x0_val_over;
%%

% Final Kernal and auxiliary vector
ST = C'*P*C; 
vT = C'*P*r1(end);


% Solve the ODEs backwards in time
[tsol1,SV1] = ode45(@(t,SV) RiccatiV2(t,SV,A,B,C,Q,R,flip(r1),tr),[0 T],[ST(1,1) ST(1,2) ST(2,2) ST(1,3) ST(1,4) ST(2,3) ST(2,4) ST(3,3) ST(3,4) ST(4,4) vT(1) vT(2) vT(3) vT(4)]);


N = length(tsol1);

SVf = zeros(N,14); 
t1 = zeros(N,1);
for k = 1:N
    SVf(k,:) = SV1(N-k+1,:); %S and v forward in time by flipping index  
    
    S1= [SV1(N-k+1,1),SV1(N-k+1,2),SV1(N-k+1,4),SV1(N-k+1,5);
              SV1(N-k+1,2),SV1(N-k+1,3),SV1(N-k+1,6),SV1(N-k+1,7);
               SV1(N-k+1,4),SV1(N-k+1,6),SV1(N-k+1,8),SV1(N-k+1,9);
               SV1(N-k+1,5),SV1(N-k+1,7),SV1(N-k+1,9),SV1(N-k+1,10)]; 
           
    K(k,:) = inv(R)*B'*S1; % calculate K forward in time
    
    t1(k) = T-tsol1(N-k+1); %time information
end


%%  Plot S, V, K

figure(1)
plot(tsol1,SVf(:,1:10))
grid on
xlabel('Time (sec)')
ylabel('S(t)')
legend('s_1','s_2','s_3','s_4','s_5','s_6','s_7','s_8','s_9','s_{10}','location','northeastoutside')
title('Kernal matrix trajectory S(t)')


figure(2)
plot(tsol1,SVf(:,11:14))
grid on
xlabel('Time (sec)')
ylabel('v(t)')
legend('v1(t)','v2(t)','v3(t)','v4(t)','location','best')
title('Modified state transition matrix trajectory v(t)')


figure(3)
plot(tsol1,K)
grid on
xlabel('Time (sec)')
ylabel('K(t)')
legend('K1(t)','K2(t)','K3(t)','K4(t)','location','best')
title('Feedback Gain K(t)')


%% Solve the plant forward in time

[tcl,xf] = ode45(@(t,x1) plant1(t,x1,A,B,R,SVf,t1),[0 T],[x0(1) x0(2) x0(3) x0(4)]); %integration forward in time

%%
s1 = interp1(t1,SVf(:,1),tcl);
s2 = interp1(t1,SVf(:,2),tcl);
s3 = interp1(t1,SVf(:,3),tcl);
s4=interp1(t1,SVf(:,4),tcl);
s5=interp1(t1,SVf(:,5),tcl);
s6=interp1(t1,SVf(:,6),tcl);
s7=interp1(t1,SVf(:,7),tcl);
s8=interp1(t1,SVf(:,8),tcl);
s9=interp1(t1,SVf(:,9),tcl);
s10=interp1(t1,SVf(:,10),tcl);
V1 = interp1(t1,SVf(:,11),tcl);
V2 = interp1(t1,SVf(:,12),tcl);
V3=interp1(t1,SVf(:,13),tcl);
V4=interp1(t1,SVf(:,14),tcl);
r = interp1(tr,r1,tcl);
N1=length(tcl);
for k=1:N1-1
    
S0 = [s1(k) s2(k) s4(k) s5(k);
      s2(k) s3(k) s6(k) s7(k);
      s4(k) s6(k) s8(k) s9(k);
      s5(k) s7(k) s9(k) s10(k)];  

K0 = inv(R)*B'*S0;

V=[V1(k);
     V2(k);
     V3(k);
     V4(k)];
 

%calculating initial value of u(t)

u0 = -(K0 - (R)^-1*B'*V*(P)^-1*V')*xf(k+1,:)' - (R)^-1*B'*V*(P)^-1*r(k);
u(k)=u0;
end
%% 
figure(4)
subplot(2,1,1)
plot(tcl,xf(:,3),'linewidth',1)
  hold on
plot(tr,r1,'linewidth',1)

xlabel('Time(s)','fontsize',11)
ylabel('X(t)','fontsize',11)
title('Closedloop State Trajectory','fontsize',11)
legend('position x','r(t)','fontsize',11, 'location','best')
grid on

subplot(2,1,2)
plot(tcl,xf(:,1:2),tcl,xf(:,4),'linewidth',1)
xlabel('Time(s)','fontsize',11)
ylabel('X(t)','fontsize',11)
legend('\theta','angular rate q','velocity vx','fontsize',11, 'location','best')
grid on
j=length(tcl)-1;
figure (5)
plot(tcl(1:j),u,'linewidth',1)
 
grid on
xlabel('Time(s)','fontsize',11)
ylabel('u(t)','fontsize',11)
title('Optimal Control','fontsize',11)

%% K_inf and feedforward gain G for simulink

K_inf=K(1,:);
G = -(C*(A-B*K_inf)^-1 *B)^-1;

