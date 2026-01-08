function AD_Bistability_m398()
%% 1. Parameters 
V1 = 0.25;
V2 = 0.11;
V3 = 2.9;      
k1 = 0.35;
k2 = 5.0;
k3 = 1.0;
k4 = 1.0;
n  = 2;         
m  = 3.98;      
epsilon = 0.1;  
sigma1  = 0.05; 
X_healthy = 0.8125;      
X_pathological = 9.0;    

%% 2. Simulation Settings 
T_max_dim = 4000;       
dt = 0.05;          
N = round(T_max_dim/dt);
time = linspace(0, T_max_dim, N+1);

drift_X = @(x,y) epsilon * (V1 + V3 * (y.^n ./ (k3^n + y.^n)) - k1 * x);
drift_Y = @(x,y) V2 + k4 * (x.^m) - k2 * y; 
X = zeros(1, N+1); 
Y = zeros(1, N+1);
X(1) = X_healthy; 
Y(1) = (V2 + k4*X(1)^m)/k2;

% Simulation Loop (Euler-Maruyama)
for i = 1:N
    fx = drift_X(X(i), Y(i));
    fy = drift_Y(X(i), Y(i));
    X(i+1) = X(i) + fx*dt + sigma1*sqrt(epsilon)*X(i)*randn*sqrt(dt);
    Y(i+1) = Y(i) + fy*dt; 
    if X(i+1) < 0, X(i+1) = 0; end
    if Y(i+1) < 0, Y(i+1) = 0; end
end

%% 4. Plotting
figure(11); clf; 
set(gcf, 'Color', 'w'); 
hold on;
h_lower = yline(X_healthy, '--', 'Color', [0.00, 0.45, 0.74], 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Healthy state (X = %.4f)', X_healthy)); 
h_upper = yline(X_pathological, '--', 'Color', [0.85, 0.33, 0.10], 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Pathological state (X = %.4f)', X_pathological));
h_path = plot(time, X, 'Color', [0.47, 0.67, 0.19], 'LineWidth', 0.8, ...
    'DisplayName', 'Sample path');
xlabel('Dimensionless Time (t)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('A\beta Concentration (X)', 'FontSize', 14, 'FontName', 'Times New Roman', 'Rotation', 0, 'HorizontalAlignment', 'right');
title('Bistability with m=3.98', 'FontSize', 14, 'FontName', 'Times New Roman');
xlim([0 T_max_dim]);
ylim([0 10]);
xtick_positions = 0:500:T_max_dim;
set(gca, 'XTick', xtick_positions);
set(gca, 'Box', 'on', 'FontSize', 12, 'LineWidth', 1.0, 'FontName', 'Times New Roman');
set(gca, 'YTick', 0:2:10);
lgd = legend([h_lower, h_path, h_upper], ...
    {'Healthy state', 'Sample path', 'Pathological state'}, ...
    'Location', 'northwest', 'EdgeColor', 'none');
set(lgd, 'FontSize', 10, 'FontName', 'Times New Roman');

hold off;

%% 5. Add text annotation with simulation info
annotation('textbox', [0.65 0.15 0.28 0.12], 'String', ...
    sprintf('Parameters (Shaheen m=3.98):\nε = %.1f, σ₁ = %.2f, m = %.2f\nX_healthy = %.4f, X_path = %.2f', ...
    epsilon, sigma1, m, X_healthy, X_pathological), ...
    'FontSize', 9, 'FitBoxToText', 'on', 'BackgroundColor', [0.95 0.95 0.95], 'FontName', 'Times New Roman');

end
