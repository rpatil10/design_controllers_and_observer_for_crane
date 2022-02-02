syms l1 l2 m1 m2 g M
m1 = 100;
m2 = 100;
M = 1000;
l1 = 20;
l2 = 10;
g = 9.81;
q0 = [2 0 deg2rad(17) 0 deg2rad(30) 0];
tspan = 0:0.1:100;
A = [0 1 0 0 0 0; 0 0 -m1*g/M 0 -m2*g/M 0; 0 0 0 1 0 0; 0 0 -((M*g)+(m1*g))/(M*l1) 0 -g*m2/(M*l1) 0; 0 0 0 0 0 1; 0 0 -m1*g/(M*l2) 0 -((M*g)+(m2*g))/(M*l2) 0];
B = [0; 1/M; 0; 1/(l1*M); 0; 1/(l2*M)];
C1 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
C2 = [0 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
C3 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 1 0];
C4 = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
D = [0; 0; 0];
Observer1 = rank([C1' A'*C1' ((A')^2)*C1' ((A')^3)*C1' ((A')^4)*C1' ((A')^5)*C1']);
Observer2 = rank([C2' A'*C2' ((A')^2)*C2' ((A')^3)*C2' ((A')^4)*C2' ((A')^5)*C2']);
Observer3 = rank([C3' A'*C3' ((A')^2)*C3' ((A')^3)*C3' ((A')^4)*C3' ((A')^5)*C3']);
Observer4 = rank([C4' A'*C4' ((A')^2)*C4' ((A')^3)*C4' ((A')^4)*C4' ((A')^5)*C4']);
SYS1 = ss(A,B,C1,D);
SYS3 = ss(A,B,C3,D);
SYS4 = ss(A,B,C4,D);
Bd = 0.1*eye(6);                
Vn = 0.01;                      
[Lue_1,P,E] = lqe(A,Bd,C1,Bd,Vn*eye(3));
[Lue_3,P,E] = lqe(A,Bd,C3,Bd,Vn*eye(3));
[Lue_4,P,E] = lqe(A,Bd,C4,Bd,Vn*eye(3));
AC1 = A-(Lue_1*C1);
AC3 = A-(Lue_3*C3);
AC4 = A-(Lue_4*C4);
sys1 = ss(AC1,[B Lue_1],C1,0);
sys3 = ss(AC3,[B Lue_3],C3,0);
sys4 = ss(AC4,[B Lue_4],C4,0);
unitStep = 0*tspan;
unitStep(200:length(tspan)) = 1;
[y1,t] = lsim(SYS1,unitStep,tspan);
[x1,t] = lsim(sys1,[unitStep;y1'],tspan);
[y3,t] = lsim(SYS3,unitStep,tspan);
[x3,t] = lsim(sys3,[unitStep;y3'],tspan);
[y4,t] = lsim(SYS4,unitStep,tspan);
[x4,t] = lsim(sys4,[unitStep;y4'],tspan);
figure();
hold on
plot(t,y1(:,1),'r','Linewidth',2)
plot(t,x1(:,1),'k--','Linewidth',1)
ylabel('State Variables')
xlabel('time(sec)')
legend('x(t)','Estimated x(t)')
title('Response for output vector at step input: (x(t)')
hold off
figure();
hold on
plot(t,y3(:,1),'r','Linewidth',2)
plot(t,y3(:,3),'b','Linewidth',2)
plot(t,x3(:,1),'k--','Linewidth',1)
plot(t,x3(:,3),'m--','Linewidth',1)
ylabel('State Variables')
xlabel('time(sec)')
legend('x(t)','theta_2(t)','Estimated x(t)','Estimated theta_2(t)')
title('Response for output vector at step input: (x(t),theta_2(t))')
hold off
figure();
hold on
plot(t,y4(:,1),'r','Linewidth',2)
plot(t,y4(:,2),'g','Linewidth',2)
plot(t,y4(:,3),'b','Linewidth',2)
plot(t,x4(:,1),'k--','Linewidth',1)
plot(t,x4(:,2),'r--','Linewidth',1)
plot(t,x4(:,3),'m--','Linewidth',1)
ylabel('State Variables')
xlabel('time(sec)')
legend('x(t)','theta_1(t)','theta_2(t)','Estimated x(t)','Estimated theta_1(t)','Estimated theta_2(t)')
title('Response for output vector at step input: (x(t),theta_1(t),theta_2(t))')
hold off
%% Linear system response
[t,q1] = ode45(@(t,q)linearObs1(t,q,Lue_1),tspan,q0);
figure();
hold on
plot(t,q1(:,1))
ylabel('state variables')
xlabel('time (sec)')
title('Linear system Observer for output vector: x(t)')
legend('x')
hold off
[t,q3] = ode45(@(t,q)linearObs3(t,q,Lue_3),tspan,q0);
figure();
hold on
plot(t,q3(:,1))
plot(t,q3(:,5))
ylabel('state variables')
xlabel('time (sec)')
title('Linear system Observer for output vector: (x(t),theta_2(t))')
legend('x','theta_2')
hold off
[t,q4] = ode45(@(t,q)linearObs4(t,q,Lue_4),tspan,q0);
figure();
hold on
plot(t,q4(:,1))
plot(t,q4(:,3))
plot(t,q4(:,5))
ylabel('state variables')
xlabel('time (sec)')
title('Linear system Observer for output vector: (x(t),theta_1(t),theta_2(t))')
legend('x','theta_1','theta_2')
hold off
%% Non-linear system observer response
[t,q1] = ode45(@(t,q)nonLinearObs1(t,q,1,Lue_1),tspan,q0);
figure();
hold on
plot(t,q1(:,1))
ylabel('state variables')
xlabel('time (sec)')
title('Non-Linear System Observer for output vector: x(t)')
legend('x')
hold off
[t,q3] = ode45(@(t,q)nonLinearObs3(t,q,1,Lue_3),tspan,q0);
figure();
hold on
plot(t,q3(:,1))
plot(t,q3(:,5))
ylabel('state variables')
xlabel('time (sec)')
title('Non-Linear System Observer for output vector: (x(t),theta_2(t))')
legend('x','theta_2')
hold off
[t,q4] = ode45(@(t,q)nonLinearObs4(t,q,1,Lue_4),tspan,q0);
figure();
hold on
plot(t,q4(:,1))
plot(t,q4(:,3))
plot(t,q4(:,5))
ylabel('state variables')
xlabel('time (sec)')
title('Non-Linear System Observer for output vector: (x(t),theta_1(t),theta_2(t))')
legend('x','theta_1','theta_2')
hold off