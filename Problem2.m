%% MAE 263F HW1 P2

%% General Spheres Falling in Viscous Fluid

% Number of nodes
N=27;
ndof = N * 2;

dt = 0.01;

RodLength = 1; % 1 meter
deltaL = RodLength / (N-1);

% Radii of spheres in meters
R = zeros(N,1); % vector of size N
R(:) = deltaL/10;
midNode = (N+1)/2;
R(midNode) = 0.025;

% Density
rho_metal = 7000;
rho_f = 1000;
rho = rho_metal - rho_f;

r0 = 0.001;
Y = 1e9;
g = 9.81;

visc = 1000;
totalTime = 50;

% Utility parameter
ne = N - 1;
EI = Y * pi * r0^4 / 4;
EA = Y * pi * r0^2;

% Geometry
nodes = zeros(N, 2);
for c=1:N % loop over all nodes
    nodes(c,1) = (c-1) * deltaL;
    nodes(c,2) = 0;
end

% Mass, M
M = zeros(ndof, ndof);
for k=1:N
    M(2*k-1, 2*k-1) = 4/3*pi*R(k)^3*rho_metal; % mass for x_k
   M(2*k, 2*k) = M(2*k-1, 2*k-1); % mass for y_k
end

% Viscous damping matrix

C = zeros(ndof,ndof);
for k = 1:N
    C(2*k-1, 2*k-1) = 6 * pi * R(k);
    C(2*k, 2*k) = C(2*k-1, 2*k-1);
end

% Weight vector, W

W = zeros(ndof, 1);
for k = 1:N
    W(2*k-1) = 0;
    W(2*k) = -4/3*pi*R(k)^3*rho*g;

end

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

for c = 2:Nsteps
    fprintf("Time = %f\n", (c-1) * dt);

    % guess
    q = q0;

    % Newton Raphson
    err = 10 * tol;
    while err > tol
        f = M/dt * ((q-q0) / dt - u0);
        J = M / dt^2;


        % linear springs
        for k=1:N-1
            xk = q(2*k-1);
            yk = q(2*k);
            xkp1 = q(2*k+1);
            ykp1 = q(2*k+2);
            l_k = deltaL;
            dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA);
            dJ = hessEs(xk, yk, xkp1, ykp1, l_k, EA);
            ind = [2*k-1, 2*k, 2*k+1, 2*k+2];
            f(ind) = f(ind) + dF;
            J(ind, ind) = J(ind, ind) + dJ;
        end

        % Bending springs
        for k=2:N-1
            xkm1 = q(2*k-3);
            ykm1 = q(2*k-2);
            xk = q(2*k-1);
            yk = q(2*k);
            xkp1 = q(2*k+1);
            ykp1 = q(2*k+2);
            curvature0 = 0;
            l_k = deltaL;
            dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
            dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
            ind = [2*k-3, 2*k-2, 2*k-1, 2*k, 2*k+1, 2*k+2];
            f(ind) = f(ind) + dF;
            J(ind, ind) = J(ind, ind) + dJ;
        end

        % Viscous force
        f = f + C * (q-q0) / dt;
        J = J + C /dt;

        % Weight
        f = f - W;

        % Update
        q = q - J \f;

        err = sum (abs(f));

    end

    % New velocity
    u = (q- q0) / dt;
    
    % Store position info
    all_mid_y(c) = q(4);

    % Store velocity info
    all_mid_v(c) = u(2*midNode);
%{
    % Plot
    figure(1);
    plot( q(1:2:end), q (2:2:end), "ro-");
    axis equal
    xlabel("c [meters]");
    ylabel("y [meters]");
    drawnow
%}

    % Update
    q0 = q;
    u0 = u;
end


% plot middle node position

figure(2)
timeArray = (1:Nsteps) * dt;
plot(timeArray, all_mid_y, "r-");
xlabel("Time, t [seconds]");
ylabel("Vertical position of middle node, y [m]");
title("Middle Node Position")

% plot middle node downward velocity
figure(3)
timeArray = (1:Nsteps) * dt;
plot(timeArray, all_mid_v, "r-");
xlabel("Time, t [seconds]");
ylabel("Vertical velocity of middle node, v [m/s]");
title("Middle Node Velocity")
