function AD_Simulations()
    AD_bistable();
    AD_Bifurcation();
    AD_SDE_Simulation();
end

%% --------------------------------------
%% PART 1: AD_bistable.m
%% --------------------------------------
% Function to visualize the phase plane and the sigmoidal behavior of amyloid
function AD_bistable()
%% General parameters (fixed for all the simulations)
V1 = 0.25;
V2 = 0.11;
k1 = 0.35;
k2 = 5.0;
k3 = 1.0;
k4 = 1.0;
n = 2; 

% Settings
opts = odeset('RelTol',1e-8);

%% FIRST SIMULATION: V3 = 2.9 (bistability)
% Variable paramters
epsilon_phase = 0.01; 
V3_bistable = 2.9; 

pars_bistable = {V1, V2, V3_bistable, k1, k2, k3, k4, epsilon_phase, n};
RHS_bistable = @(t,y) rhs(t,y,pars_bistable);

% Initial conditions for the trajectories
ICs = [...
    1.0, 0.5;    
    2.85, 0.45;  
    4.5, 0.9;    
    2.09, 0.9];  

T_phase = 2000;
x_all = cell(size(ICs,1),1);     
y_all = cell(size(ICs,1),1); 

% Computation of the trajectories
for j = 1:size(ICs, 1)
    [t, Y] = ode15s(RHS_bistable, [0 T_phase], ICs(j,:)', opts);
    x_all{j} = Y(:,1); 
    y_all{j} = Y(:,2); 
end

% Nullclines
y_vec = 0:0.01:2.5; 
x_null_bistable = (V1 + V3_bistable .* (y_vec.^n) ./ (k3^n + y_vec.^n)) ./ k1;
x_req = (k2 .* y_vec - V2) ./ k4;

%% SECOND SIMULATION: V3 = 3.05 (Fig 1 - sigmoidal)
% Variable paramters
V3_disease = 3.05; 
epsilon_fig1 = 0.01; 

pars_disease = {V1, V2, V3_disease, k1, k2, k3, k4, epsilon_fig1, n};
RHS_disease = @(t,y) rhs(t,y,pars_disease);

IC_fig1 = [2.09, 0.9];
[t_sig, Y_sig] = ode15s(RHS_disease, [0 12000], IC_fig1', opts);
x_sigmoid = Y_sig(:,1);

%% PLOTTING
% Exact computation of equilibria 
a_coeff = -k1 * k4^2;
b_coeff = -2*V2*k4*k1 + k4^2*(V1 + V3_bistable);
c_coeff = -(V2^2 + k2^2*k3^2)*k1 + 2*V2*k4*(V1 + V3_bistable);
d_coeff = V1*k2^2*k3^2 + V2^2*(V1 + V3_bistable);

x_roots = roots([a_coeff, b_coeff, c_coeff, d_coeff]);
x_eq = sort(real(x_roots(imag(x_roots)==0))); 
y_eq = (k4 .* x_eq + V2) ./ k2; 

% FIGURE 1: Phase Planes
figure(1); clf; 
set(gcf, 'Name', 'AD Bistable - Phase Planes', 'NumberTitle', 'off');

subplot(1,2,1); hold on;
plot(x_null_bistable, y_vec, 'r-', 'LineWidth', 1);      
plot(x_req, y_vec, 'b-', 'LineWidth', 1);       
colors = lines(size(ICs,1));
for j = 1:size(ICs,1)
    plot(x_all{j}, y_all{j}, 'Color', colors(j,:), 'LineWidth', 1.5);
    plot(x_all{j}(1), y_all{j}(1), 'Color', colors(j,:), 'MarkerFaceColor', colors(j,:));
end

% Stability Markers
plot(x_eq(1), y_eq(1), 'ko', 'MarkerFaceColor', 'g', 'MarkerSize', 3); % Healthy
plot(x_eq(3), y_eq(3), 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 3); % Pathological
plot(x_eq(2), y_eq(2), 'kx', 'LineWidth', 2, 'MarkerSize', 3); % Saddle
xlabel('x'); ylabel('y');
title('Phase plane');
axis([0 8 0 2]); grid on;
legend({'Nullcline x', 'Nullcline y', 'Trajectories'}, 'Location', 'best');
hold off;

subplot(1,2,2); hold on;
plot(x_null_bistable, y_vec, 'r-', 'LineWidth', 1);      
plot(x_req, y_vec, 'b-', 'LineWidth', 1);       
for j = 1:size(ICs,1)
    plot(x_all{j}, y_all{j}, 'Color', colors(j,:), 'LineWidth', 1.5);
end
plot(x_eq(1), y_eq(1), 'ko', 'MarkerFaceColor', 'g', 'MarkerSize', 3); 
plot(x_eq(3), y_eq(3), 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 3); 
plot(x_eq(2), y_eq(2), 'kx', 'LineWidth', 2, 'MarkerSize', 3); 
xlabel('x'); ylabel('y');
title('Phase plane zoom');
grid on; grid minor;
xlim([1.2 5.0]); ylim([0.25 1.0]); 
hold off;

% FIGURE 2: Time Course
figure(2); clf; hold on;
set(gcf, 'Name', 'AD Bistable - Time Course', 'NumberTitle', 'off');
plot(t_sig, x_sigmoid, 'b-', 'LineWidth', 3);

xlabel('Time (t)');
ylabel('X concentration');
title(['Sigmoidal Transition (V_3=' num2str(V3_disease) ')']);
grid on;
ylim([1 6]);
hold off;

end

function dY = rhs(t,Y,pars)
[V1, V2, V3, k1, k2, k3, k4, epsilon, n] = pars{:}; 
x = Y(1); 
y = Y(2); 

dx = epsilon * (V1 + V3 * (y^n / (k3^n + y^n)) - k1 * x);
dy = V2 + k4 * x - k2 * y;

dY = [dx; dy];
end


%% --------------------------------------
%% PART 2: AD_Bifurcation.m
%% --------------------------------------
% Function to generate the bifurcation diagram
function AD_Bifurcation()
%% Fixed parameters
V1 = 0.25;
V2 = 0.11;
k2 = 5.0;
k3 = 1.0;
k4 = 1.0;

V3_range = 0:0.002:5; 

figure(3); clf; hold on
set(gcf, 'Color', 'w', 'Name', 'AD Bifurcation Analysis', 'NumberTitle', 'off'); 

% PLOT A: Bistable (k1 = 0.35)
subplot(1, 2, 1); hold on;
k1_a = 0.35; 

[v3_s, x_s, v3_u, x_u] = calculate_branches(V1, V2, k1_a, k2, k3, k4, V3_range);

plot(v3_s, x_s, 'b.', 'MarkerSize', 6); 
if ~isempty(v3_u)
    plot(v3_u, x_u, 'g.', 'MarkerSize', 6);
end

title(['(a) k_1 = ' num2str(k1_a)], 'FontSize', 14);
xlabel('V_3', 'FontSize', 12);
ylabel('x at equilibrium', 'FontSize', 12);
grid on; grid minor;
xlim([2.0 3.8]); 
ylim([0 8]);

% PLOT B: Unimodal (k1 = 0.293)
subplot(1, 2, 2); hold on;
k1_b = 0.293; 

[v3_s, x_s, v3_u, x_u] = calculate_branches(V1, V2, k1_b, k2, k3, k4, V3_range);

plot(v3_s, x_s, 'b.', 'MarkerSize', 5); 

title(['(b) k_1 = ' num2str(k1_b)], 'FontSize', 14);
xlabel('V_3', 'FontSize', 12);
ylabel('x at equilibrium', 'FontSize', 12);
grid on;
xlim([0 5]); 
ylim([0 15]);

hold off;
end

function [v3_stable, x_stable, v3_unstable, x_unstable] = calculate_branches(V1, V2, k1, k2, k3, k4, V3_range)
    v3_stable = []; x_stable = [];
    v3_unstable = []; x_unstable = [];
    
    for V3 = V3_range
        A = -k1 * k4^2;
        B = -2*V2*k4*k1 + k4^2*(V1 + V3);
        C = -(V2^2 + k2^2*k3^2)*k1 + 2*V2*k4*(V1 + V3);
        D = V1*k2^2*k3^2 + V2^2*(V1 + V3);
        
        all_roots = roots([A, B, C, D]);
        
        real_roots = all_roots(imag(all_roots)==0 & real(all_roots)>=0);
        real_roots = sort(real_roots); 
        
        if length(real_roots) == 1
            v3_stable = [v3_stable, V3];
            x_stable = [x_stable, real_roots(1)];
        elseif length(real_roots) == 3
            v3_stable = [v3_stable, V3];
            x_stable = [x_stable, real_roots(1)];
            
            v3_unstable = [v3_unstable, V3];
            x_unstable = [x_unstable, real_roots(2)];
            
            v3_stable = [v3_stable, V3];
            x_stable = [x_stable, real_roots(3)];
        end
    end
end

%% --------------------------------------
%% PART 3: AD_SDE_Simulation.m
%% --------------------------------------
function AD_SDE_Simulation()
% Fixed Parameters 
V1 = 0.25;
V2 = 0.11;
V3 = 2.9; 
k1 = 0.35;
k2 = 5.0;
k3 = 1.0;
k4 = 1.0;
n  = 2;         
epsilon = 0.2;

% Noise intensity
sigma1 = 0.05;  
sigma2 = 0.0;   

% SDE Simulation
T_max = 4000;       
dt = 0.05;          
N = round(T_max/dt);
time = linspace(0, T_max, N+1);

X = zeros(1, N+1); 
Y = zeros(1, N+1); 

% Initial Condition: healthy state
X(1) = 1.45; 
Y(1) = (V2 + k4*X(1))/k2; 

% Simulation Loop
for i = 1:N
    x_curr = X(i);
    y_curr = Y(i);
    
    % Deterministic drift
    f_x = epsilon * (V1 + V3 * (y_curr^n / (k3^n + y_curr^n)) - k1 * x_curr);
    f_y = V2 + k4 * x_curr - k2 * y_curr;
    
    % Stochastic diffusion
    g_x = sigma1 * sqrt(epsilon) * x_curr;
    g_y = sigma2 * y_curr;
    
    % Update
    X(i+1) = x_curr + f_x * dt + g_x * (randn * sqrt(dt));
    Y(i+1) = y_curr + f_y * dt + g_y * (randn * sqrt(dt));
    
    % Positivity constraint
    if X(i+1) < 0, X(i+1) = 0; end
    if Y(i+1) < 0, Y(i+1) = 0; end
end

% Plotting
% SDE Path
figure(4); clf; 
set(gcf, 'Color', 'w', 'Name', 'AD SDE Simulation', 'NumberTitle', 'off'); 
hold on;

% 1. Plot the Sample Path (Green)
h_path = plot(time, X, 'Color', [0.47, 0.67, 0.19], 'LineWidth', 0.8); % Greenish

% 2. Plot Steady States (Reference Lines)
h_upper = yline(4.7, 'Color', [0.85, 0.33, 0.10], 'LineWidth', 1.5); 
h_lower = yline(1.45, 'Color', [0.00, 0.45, 0.74], 'LineWidth', 1.5); 

xlabel('t', 'FontSize', 14, 'Interpreter', 'tex');
ylabel('x', 'FontSize', 14, 'Interpreter', 'tex', 'Rotation', 0, 'HorizontalAlignment', 'right');

xlim([0 T_max]);
ylim([0 7]);
title('Stochastic Differential Equation Simulation');
lgd = legend([h_upper, h_path, h_lower], ...
    {'Pathological state', 'Sample path', 'Healthy state'}, ...
    'Location', 'northwest', 'EdgeColor', 'none');

hold off;
end
