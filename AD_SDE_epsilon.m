function AD_SDE_epsilon()

V1 = 0.25;
V2 = 0.11;
V3 = 2.9;      
k1 = 0.35;
k2 = 5.0;
k3 = 1.0;
k4 = 1.0;
n  = 2;         

% Noise intensity fixed
sigma1 = 0.05;  
sigma2 = 0.0;   

T_max = 2000;       
dt = 0.05;          
N = round(T_max/dt);
time = linspace(0, T_max, N+1);

% 3 scenari per Epsilon
eps_values = [0.01, 0.5, 1.0];
labels = {'\epsilon = 0.01 (paper)', '\epsilon = 0.5 (intermediate)', '\epsilon = 1.0 (no separation)'};
colors = [0 0 1; 0.85 0.33 0.10; 0.47 0.67 0.19];

X_results = zeros(3, N+1);

% Stesso rumore per tutte e 3 le simulazioni
rng(42); % seed 
dW1_common = randn(1, N) * sqrt(dt);
dW2_common = randn(1, N) * sqrt(dt);

for k = 1:length(eps_values)
    epsilon = eps_values(k);
    
    X = zeros(1, N+1); 
    Y = zeros(1, N+1); 
    
    X(1) = 1.45; 
    Y(1) = (V2 + k4*X(1))/k2; 
    
    for i = 1:N
        x_curr = X(i);
        y_curr = Y(i);
        
        f_x = epsilon * (V1 + V3 * (y_curr^n / (k3^n + y_curr^n)) - k1 * x_curr);
        f_y = V2 + k4 * x_curr - k2 * y_curr;
        
     
        g_x = sigma1 * sqrt(epsilon) * x_curr;
        g_y = sigma2 * y_curr;
        

        X(i+1) = x_curr + f_x * dt + g_x * dW1_common(i);
        Y(i+1) = y_curr + f_y * dt + g_y * dW2_common(i);
        
     
        if X(i+1) < 0, X(i+1) = 0; end
        if Y(i+1) < 0, Y(i+1) = 0; end
    end
    
  
    X_results(k, :) = X;
end

figure(13); clf; 
set(gcf, 'Color', 'w', 'Position', [100 100 900 600]); 
hold on;

for k = 1:3
    plot(time, X_results(k,:), 'Color', colors(k,:), 'LineWidth', 1.2, 'DisplayName', labels{k});
end

yline(4.7, 'k--', 'Pathological threshold', 'LineWidth', 1, 'LabelHorizontalAlignment', 'left'); 
yline(1.45, 'k--', 'Healthy baseline', 'LineWidth', 1, 'LabelHorizontalAlignment', 'left'); 

title('\epsilon impact', 'FontSize', 14);
xlabel('Time (adimensional)', 'FontSize', 12);
ylabel('Amyloid concentration', 'FontSize', 12);
legend('Location', 'best');
grid on;
ylim([0 6]);
xlim([0 T_max]);

hold off;
end