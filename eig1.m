m1 = 100;
m2 = 100;
M = 1000;
L1 = 20;
L2 = 10;
g = 9.81;
x0 = [0; 0; deg2rad(5); 0; deg2rad(5); 0]
t = 0:0.1:80;
A = [0 1 0 0 0 0; 0 0 -g*m1/M 0 -g*m2/M 0; 0 0 0 1 0 0; 0 0 -((M*g)+(m1*g))/(M*L1) 0 -g*m2/(M*L1) 0; 0 0 0 0 0 1; 0 0 -m1*g/(M*L2) 0 -((M*g)+(m2*g))/(M*L2) 0];
B = [0; 1/M; 0; 1/(L1*M); 0; 1/(L2*M)];
C = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
D = [1;0;0];
states = {'x' 'x_dot' 'theta1' 'theta1_dot' 'theta2' 'theta2_dot'};
inputs = {'u'};
outputs = {'x'; 'theta1'; 'theta2'};
sys_ss = ss(A,B,C,D,'statename',states,'inputname',inputs,'outputname',outputs);
poles = eig(A)
Q = [1/12 0 0 0 0 0; 0 8.33 0 0 0 0; 0 0 312.5 0 0 0; 0 0 0 833.33 0 0; 0 0 0 0 312.5 0; 0 0 0 0 0 833.33];
R = 0.00001;
K = lqr(A,B,Q,R)
states = {'x' 'x_dot' 'theta1' 'theta1_dot' 'theta2' 'theta2_dot'};
inputs = {'r'};
outputs = {'x'; 'theta1'; 'theta2'};
sys_ss = ss((A-B*K),B,C,D,'statename',states,'inputname',inputs,'outputname',outputs);
initial(sys_ss,x0,t)
eig(A-B*K)