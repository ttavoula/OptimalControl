clc
clear
close all
load('antx_2022-12-07-12-55-17.mat')
%%
x.value=x.value-0.02;
x.timestamp=x.timestamp-5;
x0.timestamp=x0.timestamp-5;
M.timestamp=M.timestamp-5;
q.timestamp=q.timestamp-5;
v.timestamp=v.timestamp-5;
theta.timestamp=theta.timestamp-5;
figure(1)
plot(x0.timestamp,x0.value,'--',x.timestamp,x.value)
ylabel('Position[m]')
xlabel('Time [s]')
axis([0 15 -0.3 0.3])
figure(2)
plot(M.timestamp,M.value)
ylabel('Pitch Moment')
xlabel('Time [s]')
axis([0 15 -0.15 0.15])

figure(3)
plot(q.timestamp,q.value,theta.timestamp,theta.value,v.timestamp,v.value)
ylabel('States')
xlabel('Time [s]')
legend('pitch rate q', 'angle \theta','velocity vx')
axis([0 15 -3 3])


