function AD_Bifurcation_parametersensitivity()
    clc; clear; close all;

    V1 = 0.25;
    V2 = 0.11;
    % V3 e k1 parametri variati
    k2 = 5.0;
    k3 = 1.0;
    k4 = 1.0;
    
    % Range per V3
    V3_range = 0:0.01:6; 
    
    figure(4); clf;
    set(gcf, 'Color', 'w', 'Position', [100, 100, 1200, 500]); 

    % PLOT 1: effetto clearance (k1) sulla bistabilità
    subplot(1, 2, 1); hold on;
    
    % k1 = 0.35 (standard)
    k1_base = 0.35;
    plot_bifurcation_curve(V1, V2, k1_base, k2, k3, k4, V3_range, 'b', 'Standard (k_1=0.35)');
    
    % k1 aumentato del 10% (paziente "protetto")
    % tesi: se aumentiamo la clearance proteggiamo il paziente
    k1_better = k1_base * 1.10; 
    plot_bifurcation_curve(V1, V2, k1_better, k2, k3, k4, V3_range, 'g', '+10% Clearance (k_1=0.385)');
    
    % k1 ridotto del 10% (paziente a rischio)
    k1_worse = k1_base * 0.90;
    plot_bifurcation_curve(V1, V2, k1_worse, k2, k3, k4, V3_range, 'r', '-10% Clearance (k_1=0.315)');
    
    title('Clearance (k_1) sensibility', 'FontSize', 14);
    xlabel('Feedback Calcium', 'FontSize', 12);
    ylabel('Amyloid equilibrium', 'FontSize', 12);
    legend('Location', 'NorthWest');
    grid on; 
    xlim([1.5 4.5]); ylim([0 10]);
    
    % PLOT 2: perdita della bistabilità
    subplot(1, 2, 2); hold on;
    
    plot_bifurcation_curve(V1, V2, k1_base, k2, k3, k4, V3_range, 'b', 'Bistable');

    k1_critical = 0.22; % Valore ipotetico di perdita bistabilità, possiamo provarne altri (paper k1 = 0.293)
    plot_bifurcation_curve(V1, V2, k1_critical, k2, k3, k4, V3_range, 'm', ['Monostable (k_1=' num2str(k1_critical) ')']);
    
    title('Bistability collapse', 'FontSize', 14);
    xlabel('Feedback Calcium', 'FontSize', 12);
    ylabel('Amyloid equilibrium', 'FontSize', 12);
    legend('Location', 'NorthWest');
    grid on;
    xlim([0 5]); ylim([0 15]);
end


function plot_bifurcation_curve(V1, V2, k1, k2, k3, k4, V3_range, color, label_text)
    v3_stable = []; x_stable = [];
    v3_unstable = []; x_unstable = [];
    
    for V3 = V3_range
        
        target_fun = @(x) V1 + V3 * (((V2+x)/k2).^2 ./ (k3^2 + ((V2+x)/k2).^2)) - k1*x;
        
        x_vals = 0:0.1:20;
        f_vals = arrayfun(target_fun, x_vals);
        
  
        sign_changes = find(diff(sign(f_vals)) ~= 0);
        
        current_roots = [];
        for idx = sign_changes
 
            root = fzero(target_fun, [x_vals(idx), x_vals(idx+1)]);
            current_roots = [current_roots, root];
        end
        

        if length(current_roots) >= 1
          
            v3_stable = [v3_stable, V3];
            x_stable = [x_stable, min(current_roots)];
            
            if length(current_roots) >= 3
             
                roots_sorted = sort(current_roots);
                
                v3_unstable = [v3_unstable, V3];
                x_unstable = [x_unstable, roots_sorted(2)];
                
                v3_stable = [v3_stable, V3];
                x_stable = [x_stable, roots_sorted(3)];
            elseif length(current_roots) == 2
                 
                 v3_stable = [v3_stable, V3];
                 x_stable = [x_stable, max(current_roots)];
            end
        end
    end
    
  
    plot(v3_stable, x_stable, '.', 'Color', color, 'MarkerSize', 5, 'DisplayName', label_text);
    
   
    if ~isempty(v3_unstable)
        plot(v3_unstable, x_unstable, '.', 'Color', color, 'MarkerSize', 3, 'HandleVisibility', 'off');
    end
end