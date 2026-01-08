function AD_SDE_Milstein()

V1 = 0.25;
V2 = 0.11;
V3 = 2.9;       
k1 = 0.35;
k2 = 5.0;
k3 = 1.0;
k4 = 1.0;
n  = 2;         

epsilon = 1.0;  % No separazione temporale
sigma1 = 0.8;   % Aumentato il rumore
sigma2 = 0.0;   

T_max = 100;    
dt = 0.5;       % Aumentato
N = round(T_max/dt);

time = linspace(0, T_max, N+1);

X_euler = zeros(1, N+1); 
X_milstein = zeros(1, N+1);

X_euler(1) = 1.5; 
X_milstein(1) = 1.5;

dW_array = randn(1, N) * sqrt(dt);

for i = 1:N
 
    x_eu = X_euler(i);
    x_mil = X_milstein(i);
    
    y_curr = 0.5; 
    Hill = y_curr^n / (k3^n + y_curr^n);
    
    % EULER
    drift_eu = epsilon * (V1 + V3 * Hill - k1 * x_eu);
    g_eu = sigma1 * sqrt(epsilon) * x_eu;
    
    X_euler(i+1) = x_eu + drift_eu*dt + g_eu*dW_array(i);
    
    % MILSTEIN
    drift_mil = epsilon * (V1 + V3 * Hill - k1 * x_mil);
    g_mil = sigma1 * sqrt(epsilon) * x_mil;
    g_prime = sigma1 * sqrt(epsilon);
    
    % Correzione 
    correction = 0.5 * g_mil * g_prime * (dW_array(i)^2 - dt);
    
    X_milstein(i+1) = x_mil + drift_mil*dt + g_mil*dW_array(i) + correction;
    
    % Opzionale check 
    if X_euler(i+1) > 20, X_euler(i+1) = 20; end
    if X_milstein(i+1) > 20, X_milstein(i+1) = 20; end
end

figure(12); clf; set(gcf, 'Color', 'w');

subplot(2,1,1);
plot(time, X_euler, 'b.-', 'LineWidth', 1, 'MarkerSize', 10); hold on;
plot(time, X_milstein, 'r.--', 'LineWidth', 1.5, 'MarkerSize', 10);
legend('Euler (raw)', 'Milstein (corretto)', 'Location', 'Best');
title(['Comparison with time step dt = ' num2str(dt)]);
ylabel('A\beta concentration'); grid on;

subplot(2,1,2);
area(time, abs(X_euler - X_milstein), 'FaceColor', 'k', 'FaceAlpha', 0.2);
title('Error');
ylabel('Difference'); xlabel('Time'); grid on;
end