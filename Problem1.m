%% 263F_HW1_P1
%% Jake Kremer
%% Fall 2024

%% 3 Falling Spheres in Viscous Fluid

% Select implicit (1) or explicit (2)
chooseMethod = 1; % choose method

% Number of nodes
N=3;
ndof = N * 2;

dt = 0.01;

RodLength = 1; % 1 meter
deltaL = RodLength / (N-1);

% Radii of spheres in meters
R1 = 0.005;
R2 = 0.025;
R3 = 0.005;

% Density
rho_metal = 7000;
rho_f = 1000;
rho = rho_metal - rho_f;

r0 = 0.001;
Y = 1e9;
g = 9.81;

visc = 1000;
totalTime = 10;

% Utility parameter
ne = N - 1;
EI = Y * pi * r0^4 / 4;
EA = Y * pi * r0^2;

% Geometry
nodes = zeros(N, 2);
for c=1:N % loop pver all nodes
    nodes(c,1) = (c-1) * deltaL;
    nodes(c,2) = 0;
end

% Mass, M
M = zeros(ndof, ndof);
M(1,1) = 4/3*pi*R1^3*rho_metal;
M(2,2) = 4/3*pi*R1^3*rho_metal;
M(3,3) = 4/3*pi*R2^3*rho_metal;
M(4,4) = 4/3*pi*R2^3*rho_metal;
M(5,5) = 4/3*pi*R3^3*rho_metal;
M(6,6) = 4/3*pi*R3^3*rho_metal;

% Viscous damping matrix

C = zeros(6,6);
C1 = 6 * pi * visc *R1;
C2 = 6 * pi * visc *R2;
C3 = 6 * pi * visc *R3;

C(1,1) = C1;
C(2,2) = C1;
C(3,3) = C2;
C(4,4) = C2;
C(5,5) = C3;
C(6,6) = C3;

% Weight vector, W

W = zeros(ndof, 1);
W(2) = -4/3*pi*R1^3*rho*g;
W(4) = -4/3*pi*R2^3*rho*g;
W(6) = -4/3*pi*R3^3*rho*g;

% initial DOF
q0 = zeros(ndof, 1);
for c=1:N % loop over nodes
    q0(2*c-1) = nodes(c,1);
    q0(2*c) = nodes(c,2);
end

u0 = zeros(ndof, 1); % initial velocity

% tolerance
tol = EI/RodLength^2 * 1e-3; % small enough force that it can be neglected

% time marching scheme
Nsteps = round(totalTime / dt);

% storage for y-position of middle node
all_mid_y = zeros(Nsteps, 1);

% storage for y-velocity of middle node
all_mid_v = zeros(Nsteps,1);

% Variables to store positions at t = 0.01, t = 0.05, t = 0.1, t = 1, and t = 10
q_t001 = [];
q_t005 = [];
q_t01 = [];
q_t1 = [];
q_t10 = [];

% Initial plot at t = 0
figure;
subplot(2,3,1); % Create a 2x3 subplot to display all the plots together
plot(q0(1:2:end), q0(2:2:end), 'ro-');
axis equal;
xlabel('x [meters]');
ylabel('y [meters]');
title('t = 0 seconds');
grid on;

if chooseMethod == 1
    for c = 2:Nsteps
        currentTime = (c-1) * dt;
        fprintf("Time = %f\n", currentTime);
    
        % Guess
        q = q0;
    
        % Newton-Raphson iteration
        err = 10 * tol;
        while err > tol
            f = M/dt * ((q-q0) / dt - u0);
            J = M / dt^2;
    
            % Linear spring between nodes 1 and 2
            xk = q(1);
            yk = q(2);
            xkp1 = q(3);
            ykp1 = q(4);
            l_k = deltaL;
            dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA);
            dJ = hessEs(xk, yk, xkp1, ykp1, l_k, EA);
            f(1:4) = f(1:4) + dF;
            J(1:4, 1:4) = J(1:4, 1:4) + dJ;
    
            % Linear spring between nodes 2 and 3
            xk = q(3);
            yk = q(4);
            xkp1 = q(5);
            ykp1 = q(6);
            l_k = deltaL;
            dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA);
            dJ = hessEs(xk, yk, xkp1, ykp1, l_k, EA);
            f(3:6) = f(3:6) + dF;
            J(3:6, 3:6) = J(3:6, 3:6) + dJ;
    
            % Bending spring at node 2
            xkm1 = q(1);
            ykm1 = q(2);
            xk = q(3);
            yk = q(4);
            xkp1 = q(5);
            ykp1 = q(6);
            curvature0 = 0;
            l_k = deltaL;
            dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
            dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
            f(1:6) = f(1:6) + dF;
            J(1:6, 1:6) = J(1:6, 1:6) + dJ;
    
            % Viscous force
            f = f + C * (q-q0) / dt;
            J = J + C /dt;
    
            % Weight
            f = f - W;
    
            % Update
            q = q - J \ f;
    
            err = sum (abs(f));
        end
    
        % New velocity
        u = (q - q0) / dt;
    
        % Store positions at t = 0.01, t = 0.05, t = 0.1, t = 1, and t = 10
        if abs(currentTime - 0.01) < 1e-6
            q_t001 = q;  % Store positions at t = 0.01
        elseif abs(currentTime - 0.05) < 1e-6
            q_t005 = q;  % Store positions at t = 0.05
        elseif abs(currentTime - 0.1) < 1e-6
            q_t01 = q;  % Store positions at t = 0.1
        elseif abs(currentTime - 1) < 1e-6
            q_t1 = q;  % Store positions at t = 1
        elseif abs(currentTime - 10) < 1e-6
            q_t10 = q;  % Store positions at t = 10
        end
    
    
        % Store position info
        all_mid_y(c) = q(4);
    
        % Store velocity info
        all_mid_v(c) = u(4);
    
        % Plot structure shape
        
        figure(1)
        plot(q(1:2:end), q(2:2:end), "ro-");
        axis equal
        xlabel("c [meters]");
        ylabel("y [meters]");
        drawnow
        
    
        % Update
        q0 = q;
        u0 = u;
    end
end

if chooseMethod == 2

    dt = 0.0001;

    for c = 2:Nsteps
    fprintf("Time = %f\n", (c-1) * dt);

    % Calculate forces (spring and damping)
    f = zeros(ndof, 1);  % Initialize force vector

    % Linear spring force between nodes 1 and 2
    xk = q0(1); yk = q0(2);
    xkp1 = q0(3); ykp1 = q0(4);
    l_k = deltaL;
    dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA);
    f(1:4) = f(1:4) + dF;

    % Linear spring force between nodes 2 and 3
    xk = q0(3); yk = q0(4);
    xkp1 = q0(5); ykp1 = q0(6);
    dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA);
    f(3:6) = f(3:6) + dF;

    % Add viscous forces
    f = f + C * u0;

    % Add weight forces
    f = f - W;

    % Velocity update (explicitly)
    u = u0 + dt * (f ./ diag(M));

    % Position update (explicitly)
    q = q0 + dt * u;

    % Store position and velocity of the middle node
    all_mid_y(c) = q(4);   % Middle node y-position
    all_mid_v(c) = u(4);   % Middle node y-velocity

    % Update the state for the next time step
    q0 = q;
    u0 = u;
    end
end

% Plot at t = 0.01
if ~isempty(q_t001)
    subplot(2,3,2);
    plot(q_t001(1:2:end), q_t001(2:2:end), 'ro-');
    axis equal;
    xlabel('x [meters]');
    ylabel('y [meters]');
    title('t = 0.01 seconds');
    grid on;
end

% Plot at t = 0.05
if ~isempty(q_t005)
    subplot(2,3,3);
    plot(q_t005(1:2:end), q_t005(2:2:end), 'ro-');
    axis equal;
    xlabel('x [meters]');
    ylabel('y [meters]');
    title('t = 0.05 seconds');
    grid on;
end

% Plot at t = 0.1
if ~isempty(q_t01)
    subplot(2,3,4);
    plot(q_t01(1:2:end), q_t01(2:2:end), 'ro-');
    axis equal;
    xlabel('x [meters]');
    ylabel('y [meters]');
    title('t = 0.1 seconds');
    grid on;
end

% Plot at t = 1
if ~isempty(q_t1)
    subplot(2,3,5);
    plot(q_t1(1:2:end), q_t1(2:2:end), 'ro-');
    axis equal;
    xlabel('x [meters]');
    ylabel('y [meters]');
    title('t = 1 seconds');
    grid on;
end

% Plot at t = 10
if ~isempty(q_t10)
    subplot(2,3,6);
    plot(q_t10(1:2:end), q_t10(2:2:end), 'ro-');
    axis equal;
    xlabel('x [meters]');
    ylabel('y [meters]');
    title('t = 10 seconds');
    grid on;
end

% plot middle node position

figure(2)
timeArray = (1:Nsteps) * dt;
plot(timeArray, all_mid_y, "b-");
xlabel("Time, t [seconds]");
ylabel("Vertical position of middle node, y [m]");
title("Middle Node Position")

% plot middle node downward velocity
figure(3)
timeArray = (1:Nsteps) * dt;
plot(timeArray, all_mid_v, "k-");
xlabel("Time, t [seconds]");
ylabel("Vertical velocity of middle node, v [m/s]");
title("Middle Node Velocity")
