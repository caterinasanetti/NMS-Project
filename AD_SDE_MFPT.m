function AD_SDE_MFPT()
    clc; clear; close all;

    V1 = 0.25;
    V2 = 0.11;
    V3 = 2.9;       
    k1 = 0.35;
    k2 = 5.0;
    k3 = 1.0;
    k4 = 1.0;
    n  = 2;         
    epsilon = 0.5;  % Ridotto scala temporale
    
    % sigma=0.05 il salto è troppo raro per vederlo in T=4000
    % Aumentato sigma per vedere le transizioni in tempi ok
    sigma1 = 0.25;  
    sigma2 = 0.05;   

    num_simulations = 500; % Numero di simulazioni
    threshold = 4.0;       % Soglia (superata barriera instabile)
    
    T_max = 5000;          % Tempo massimo di osservazione
    dt = 0.05;          
    N = round(T_max/dt);
    
    fpt_results = NaN(num_simulations, 1);
    
    
    % Monte Carlo
    for m = 1:num_simulations
        
        X = 1.45; % Stato Sano
        Y = (V2 + k4*X)/k2;
        
        passed = false; % Flag: si è ammalato?
        
        for i = 1:N
        
            Hill = Y^n / (k3^n + Y^n);
            f_x = epsilon * (V1 + V3 * Hill - k1 * X);
            f_y = V2 + k4 * X - k2 * Y;
            
            g_x = sigma1 * sqrt(epsilon) * X;
            g_y = sigma2 * Y;
            
  
            dW1 = randn * sqrt(dt);
            dW2 = randn * sqrt(dt);
            
          
            X = X + f_x * dt + g_x * dW1;
            Y = Y + f_y * dt + g_y * dW2;
            
       
            if X < 0, X = 0; end
            if Y < 0, Y = 0; end
            
            % MFPT
            % Controlla se ha superato la soglia patologica
            if X > threshold
                fpt_results(m) = i * dt; % Salva il tempo corrente
                passed = true;
                break; 
            end
        end
    
    end

    % Filtra i pazienti che NON si sono ammalati
    valid_times = fpt_results(~isnan(fpt_results));
    sick_count = length(valid_times);
    mean_fpt = mean(valid_times);
    
    fprintf('------------------------------------------------\n');
    fprintf('MFPT\n');
    fprintf('------------------------------------------------\n');
    fprintf('Total: %d\n', num_simulations);
    fprintf('Patients developing the disease within T=%d: %d (%.1f%%)\n', T_max, sick_count, (sick_count/num_simulations)*100);
    fprintf('Mean First Passage Time (MFPT): %.2f unit of time (adimensional)\n', mean_fpt);
    fprintf('------------------------------------------------\n');

    figure(14); clf; 
    set(gcf, 'Color', 'w');
    
    histogram(valid_times, 30, 'FaceColor', [0.85 0.33 0.10], 'EdgeColor', 'k', 'Normalization', 'pdf');
    hold on;
    
    xline(mean_fpt, 'k--', ['MFPT = ' num2str(mean_fpt, '%.1f')], 'LineWidth', 2);
    
    title('Distribution of first passage times', 'FontSize', 14);
    xlabel('Time to disease onset', 'FontSize', 12);
    ylabel('Probability density', 'FontSize', 12);
    grid on;
    
end