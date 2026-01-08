function AD_SDE_negconcentration()
    clc; clear; close all;

    % Extreme test settings
    V1 = 0.25;
    V2 = 0.11;
    V3 = 2.9;       
    k1 = 0.35;
    k2 = 5.0;
    k3 = 1.0;
    k4 = 1.0;
    n  = 2;         
    epsilon = 0.5;
    
    neg_violations_X = 0; 
    neg_violations_Y = 0; 
    
    sigma1 = 3.5;  % Aumentato da 0.6 a 3.5
    sigma2 = 1.0;  % Aumentato per testare anche il calcio
    
    % Aumentato passo temporale
    dt = 0.25;          
    T_max = 4000;       
    N = round(T_max/dt);
    
    time = linspace(0, T_max, N+1);
    X = zeros(1, N+1); 
    Y = zeros(1, N+1); 
    
    % Parte vicino allo zero
    X(1) = 0.1; 
    Y(1) = 0.1; 
    
    for i = 1:N
        x_curr = X(i);
        y_curr = Y(i);
        

        f_x = epsilon * (V1 + V3 * (y_curr^n / (k3^n + y_curr^n)) - k1 * x_curr);
        f_y = V2 + k4 * x_curr - k2 * y_curr;
        

        g_x = sigma1 * sqrt(epsilon) * x_curr;
        g_y = sigma2 * y_curr;
        
       
        noise_step_X = randn * sqrt(dt);
        noise_step_Y = randn * sqrt(dt);
        
        X(i+1) = x_curr + f_x * dt + g_x * noise_step_X;
        Y(i+1) = y_curr + f_y * dt + g_y * noise_step_Y;
        
        % Conteggio violazioni
        if X(i+1) < 0
            neg_violations_X = neg_violations_X + 1;
        end
            
        if Y(i+1) < 0
            neg_violations_Y = neg_violations_Y + 1;
        end
        
        % IMPORTANTE, ma ora ignorato
        %if X(i+1) < 0, X(i+1) = 1e-6; end % possible reset a un valore piccolo positivo
        %if Y(i+1) < 0, Y(i+1) = 1e-6; end
    end
    
    perc_X = (neg_violations_X / N) * 100;
    perc_Y = (neg_violations_Y / N) * 100;
    
    fprintf('------------------------------------------------\n');
    fprintf('Concentrations test \n');
    fprintf('------------------------------------------------\n');
    fprintf('Parameters: Sigma1=%.1f, dt=%.2f\n', sigma1, dt);
    fprintf('Amyloid violations (X < 0): %d (%.2f%%)\n', neg_violations_X, perc_X);
    fprintf('Calcium violations (Y < 0): %d (%.2f%%)\n', neg_violations_Y, perc_Y);
    fprintf('------------------------------------------------\n');
    

figure(10); clf; 

% Grafico Amiloide (X) 
subplot(2,1,1);
plot(time, X, 'r', 'LineWidth', 1); hold on;
area(time, min(X, 0), 'FaceColor', 'r', 'FaceAlpha', 0.3); 
yline(0, 'k-', 'LineWidth', 2);                           
title('Amyloid (X)');
ylabel('Concentration'); grid on;
ylim([min(X)*1.1, max(X)*1.1]); % adatta gli assi includendo i negativi

% Grafico Calcio (Y)
subplot(2,1,2);
plot(time, Y, 'b', 'LineWidth', 1); hold on;
area(time, min(Y, 0), 'FaceColor', 'b', 'FaceAlpha', 0.3); 
yline(0, 'k-', 'LineWidth', 2);                            
title('Calcium (Y)');
ylabel('Concentration'); xlabel('Tempo'); grid on;
ylim([min(Y)*1.1, max(Y)*1.1]); 