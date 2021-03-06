%%
syms w h r m g k c l_nat Ix Iy Iz real
syms rho cd A real % constants for drag
syms x y z alpha beta gamma      real
syms t                       real
syms q0 q1 q2 q3 q0d q1d q2d q3d real
syms omgx omgy omgz omgxd omgyd omgzd real
syms xd yd zd alphad betad gammad   real
syms alphadd betadd gammadd real
syms Talpha Tbeta Tgamma real

%% a. Find k and l_nat

l_nat = h/2;
k = 5.8*m*g/2/(h-h*l_nat/sqrt(h^2+(w/2)^2));
c = 2*sqrt(k*m)/10;

% theta_0 = atan(w/2/h);
% eqns_1 = (k*(h/cos(theta_0) - l_nat)*cos(theta_0)*2 - m*g)/m - 4.8*g == 0;
% eqns_2 = 0 + 1/2*k*(h/cos(theta_0) - l_nat)^2*2 - (h + w/2)*m*g - 1/2*k*(sqrt(2)*w/2 - l_nat)^2*2 == 0;
% eqns = [eqns_1;eqns_2];
% var = [k;l_nat];
% param = [m g w h];
% param_value = [420 9.81 18 25];
% eqns = subs(eqns,param,param_value);
% var = subs(var,param,param_value);
% [k,l_nat] = solve(eqns,var);
% %c = 2*sqrt(k*m)/10;

%% b. Build equation of motion
% Linear forces
u = [x;y;z];
ud = [xd;yd;zd];
M = diag([m m m]);
F = [0;0;-m*g];

% Euler parameters and angular velocities
omgd = [omgxd; omgyd; omgzd];
I = diag([Ix Iy Iz]);
rot_mat = [q0^2+q1^2-q2^2-q3^2, 2*(q1*q2-q0*q3)    , 2*(q1*q3-q0*q2);
           2*(q1*q2-q0*q3)    , q0^2-q1^2+q2^2-q3^2, 2*(q2*q3-q0*q1);
           2*(q1*q3-q0*q2)    , 2*(q2*q3-q0*q1)    , q0^2-q1^2-q2^2+q3^2]';

%euler = [alpha; beta; gamma];
%eulerd = [alphad; betad; gammad];
%rel_com = [-0.01,0.01,-0.1]; % center of mass relative to sphere center

com_deviate = [0;0;0];
spring1_attach = [x;y;z] + rot_mat*([0;r;0] - com_deviate);
spring2_attach = [x;y;z] + rot_mat*([0;-r;0] - com_deviate);
spring1_vec = [0;w/2;h] - spring1_attach; % spring expressed as a vector
spring2_vec = [0;-w/2;h] - spring2_attach;

% spring
C_spring1 = norm(spring1_vec) - l_nat;
C_spring2 = norm(spring2_vec) - l_nat;
C_spring1_jac = simplify(jacobian(C_spring1,u))';
C_spring2_jac = simplify(jacobian(C_spring2,u))';
F_spring1 = k*C_spring1;
F_spring2 = k*C_spring2;

% damping
spring1_vec = spring1_vec/norm(spring1_vec); % normalize to unit vector
spring2_vec = spring2_vec/norm(spring2_vec);
spring1_vel = (ud'*spring1_vec)*spring1_vec; % project velocity onto direction of spring
spring2_vel = (ud'*spring2_vec)*spring2_vec;
F_damping1 = spring1_vel*c; 
F_damping2 = spring2_vel*c;

% drag
drag = -1/2*rho*A*cd*norm(ud)*ud;

% Torque



udd = simplify(M\(F - C_spring1_jac*F_spring1 - C_spring2_jac*F_spring2 - F_damping1 - F_damping2 + drag));

%%
m = 420;
g = 9.81;
h = 25;
w = 18;
r = 1;
rho = 1.25;
A = pi;
cd = 0.5;
udd = subs(udd);

% Time step parameters
step_size = 0.01;
total_t = 5;
t = 0:step_size:total_t;

% Initialize motion
X_init = [0;0;0;0;0;0;1;0;0;0];
[X,Xd] = initializeMotion(X_init,udd,t);
X = RK4(X,udd,step_size,t);

%% Functions

function [X,Xd] = initializeMotion(X_init,qdd,t)    

    X = zeros(6,length(t));
    Xd = zeros(6,length(t));
    x = X_init(1);
    y = X_init(2);
    z = X_init(3);
    xd = X_init(4);
    yd = X_init(5);
    zd = X_init(6);
    q0 = X_init(7);
    q1 = X_init(8);
    q2 = X_init(9);
    q3 = X_init(10);

    qdd = double(subs(qdd));
    
    X(:,1) = [x; y; z; xd; yd; zd];
    Xd(:,1) = [xd; yd; zd; qdd];
end

function Xd = findGradient(X,qdd)
    x = X(1);
    y = X(2);
    z = X(3);
    xd = X(4);
    yd = X(5);
    zd = X(6);
    q0 = 0;
    q1 = 0;
    q2 = 0;
    q3 = 0;

    qdd = double(subs(qdd));
    Xd = [xd; yd; zd; qdd];
end

function X = RK4(X,qdd,step_size,t)  
    for i = 2:length(t)
        k1 = findGradient(X(:,i-1),qdd);
        k2 = findGradient(X(:,i-1) + step_size*k1/2,qdd);
        k3 = findGradient(X(:,i-1) + step_size*k2/2,qdd);
        k4 = findGradient(X(:,i-1) + step_size*k3,qdd);
        X(:,i) = X(:,i-1) + step_size*(k1 + 2*k2 + 2*k3 + k4)/6;        
    end
end

function r = rotx(phi)
    % Rotation matrix x
    r = [ ...
        1           0           0;...
        0           cos(phi)    -sin(phi);...
        0           sin(phi)    cos(phi)];
end
function r = roty(phi)
% Rotation matrix y
r = [ ...
    cos(phi)            0           sin(phi);...
    0                   1           0;...
    -sin(phi)           0           cos(phi)];
end

function r = rotz(phi)
% Rotation matrix z
r = [ ...
    cos(phi)            -sin(phi)           0;...
    sin(phi)            cos(phi)            0;...
    0                   0                   1];
end

function dist = point2line(point,line_vec)
    %line2point_vec = 



end