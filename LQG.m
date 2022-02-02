syms m1 m2 g M l1 l2 x dx 
m1 = 100;
m2 = 100;
M = 1000;
l1 = 20;
l2 = 10;
g = 9.81;
tspan = 0:0.1:100;
q_init = [2 0 deg2rad(0) 0 deg2rad(0) 0];
A = [0 1 0 0 0 0; 0 0 -m1*g/M 0 -m2*g/M 0; 0 0 0 1 0 0; 0 0 -((M*g)+(m1*g))/(M*l1) 0 -g*m2/(M*l1) 0; 0 0 0 0 0 1; 0 0 -m1*g/(M*l2) 0 -((M*g)+(m2*g))/(M*l2) 0];
B = [0; 1/M; 0; 1/(l1*M); 0; 1/(l2*M)];
C1 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
D = [1;0;0];
sys = ss(A,B,C1,D);
Q = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
R = 0.1;
[K,S,P] = lqr(A,B,Q,R);
sys = ss(A-B*K,B,C1,D);
Bd = 0.1*eye(6);                
Vn = 0.01;                      
[Lue,P,E] = lqe(A,Bd,C1,Bd,Vn*eye(3)); 
AC1 = A-(Lue*C1);
sys_ob = ss(AC1,[B Lue],C1,0);
[t,q1] = ode45(@(t,q)nonLinearObs1(t,q,-K*q,Lue),tspan,q_init);
figure();
hold on
plot(t,q1(:,1))
ylabel('state variable')
xlabel('time (sec)')
title('Non-Linear System LQG')
legend('x')
hold off