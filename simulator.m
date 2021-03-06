function simulator
close all;
clear all;
warning off all;
global area a p pf h It t count

load('tensor.mat');
%initialization
p = 48/(30*50*2); %density in g/cm3
a = 50;%unit cm
b = 30; %size
t = 2; %thickness
pf = 0.0012; %air density unit g/cm3
N = 200;%total frames
count = 0;

%mass
area = a*b;
m = 48;

%add moment of inertia of the body here
newI = diag([m*(b^2+t^2)/12 m*(a^2+t^2)/12 m*(a^2+b^2)/12]);

%add the calculated Kirchhoff tensor here
M = eye(3).*m;
Ib = newI+pf*KT(1:3,1:3); 
M = M+pf*KT(4:6,4:6);

%get I matrix
I(1:3,1:3) = Ib;
I(4:6,4:6) = M;
It = I;

% initial condition

%initial release angle
A = [0, pi/3, 0];
a1 = A(1);
a2 = A(2);
a3 = A(3);

R1 = [1, 0, 0;
    0, cos(a1), -sin(a1);
    0, sin(a1), cos(a1)];

R2 = [cos(a2), 0, sin(a2);
    0, 1, 0;
    -sin(a2), 0, cos(a2)];

R3 = [cos(a3), -sin(a3), 0;
    sin(a3), cos(a3), 0;
    0, 0, 1];

Q = zeros(N,3); %angle
EU = zeros(N,3); %euler angle
DD = cell(N,1); %trajectory data
P = zeros(N,3); %positions in trajectory

h = 0.005;%0.005 time step

tic;
srcp = zeros(6,1); %initial velocity
G = zeros(4,4); %initial state: positon and orientation
G(1:3,1:3)=R1*R2*R3; %release angle
G(1:3,4)=[0.1; 0.1;  220]; %released height
G(4,4)=1;


for k=1:N
    %integrator
    [srcp, G] = compute(srcp, G);
    
    P(k,:)=G(1:3,4)';
    Q(k,:)=srcp(4:6)';
    EU(k,:)=srcp(1:3)';
    DD{k} = G;
end
time = toc;
% save trajetory data
save('Traj.mat', 'DD');
% DrawTraj(DD);

% fprintf('[t=%f]',time);
% title('Trajectory');
% plot3(P(:,1),P(:,2),P(:,3),'linewidth',2);
% xlabel('x');
% ylabel('y');
% zlabel('z');
% axis equal;
end

function [A, G] = compute(srcp,G)

global h It

SE = se3_ad(h*srcp);
M = SE'*It*srcp + fext(srcp, G);
src = srcp + h*It\M;

options=optimset('TolFun',1e-15,'MaxFunEvals', 6e+5,'MaxIter',9e+3,'Display','off');
A = fsolve(@(src)Dyn(src, srcp, G), src, options);

G = G * se3_cay(h*A);
end

function [f] = Dyn(src,srcp,G)

global h It
S = se3_tln(h*src);
Q = se3_tln(-1.0*h*srcp);

f = S'*It*src-Q'*It*srcp-h*fext(srcp, G);
end

function [SY] = se3_cay(src)
w = src(1:3);
v = src(4:6);
B = 2*(2*eye(3) + so3_ad(w))/(4 + norm(w)*norm(w));
SY = zeros(4,4);

SY(1:3,1:3)=so3_cay(w);
SY(1:3,4) = B*v;
SY(4,4) = 1.0;
end

function [SO] = so3_cay(w)
wt = so3_ad(w);
SO = eye(3)+4.0*(wt+wt*wt/2)/(4+norm(w)*norm(w));
end

function [SC] = se3_tln(src)
SC = eye(6) - 0.5*se3_ad(src) +se3_ad(src)*se3_ad(src)./12;
end

function [SE] = se3_ad(src)
w = src(1:3);
v = src(4:6);
n = length(src);
SE = zeros(n,n);

SE(1:3, 1:3) = so3_ad(w);
SE(4:6, 1:3) = so3_ad(v);
SE(4:6, 4:6) = so3_ad(w);
end

function [SA] = so3_ad(w)
T = zeros(3);
T(1,2) = -w(3);
T(1,3) = w(2);
T(2,1) = w(3);
T(2,3) = -w(1);
T(3,1) = -w(2);
T(3,2) = w(1);
SA = T;
end

function [FG] = fext(srcp, G)
% for external force
% added yours here
global area p pf t
FG = zeros(6,1);
g = [0; 0; -980];
r = [0.0;0.0;0.1];%center of buoyancy
R = G(1:3,1:3);

%gravity
FG(1:3) = area*t*pf*cross(r', R'*g);
FG(4:6) = area*t*(p-pf)*R'*g;
% FG = FG + force(srcp, G,u);
end

