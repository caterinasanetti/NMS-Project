function AD_SSA_Simulations()
    run_AD_SSA();
    AD_DualNoise_Reproduction();
    AD_Reproduction_Fig16();
end

%% --------------------------------------
%% SSA WITH DIRECT METHOD
%% --------------------------------------
function run_AD_SSA()
% Launcher script for the Stochastic Simulation Algorithm (SSA).

% Fixed parameters
V1 = 0.25;  
V2 = 0.11;  
V3 = 2.9;   
k1 = 0.35;  
k2 = 5.0;   
k3 = 1.0;   
k4 = 1.0;   
n  = 2;     

epsilon = 0.01;
Omega = 20; 

% Reaction network definition
vMinus = [0 0; 0 0; 1 0; 0 0; 0 0; 0 1];
vPlus = [1 0; 1 0; 0 0; 0 1; 0 1; 0 0];

c = [
    epsilon * V1 * Omega,  
    0,                     
    epsilon * k1,          
    V2 * Omega,            
    k4,                    
    k2                     
];

HillParams.V3 = V3;
HillParams.K3 = k3;
HillParams.n  = n;
HillParams.epsilon = epsilon;

% Initialization
Conc_Start = [1.45, 0.35]; 
InitialState = round(Conc_Start * Omega);
tMax = 4000; 

fprintf('Starting Gillespie Simulation (Omega=%d)...\n', Omega);
[T, Dyn] = simDM_AD_Variant2(vMinus, vPlus, c, InitialState, tMax, Omega, HillParams);

Dyn_Conc = Dyn / Omega; 

% FIGURE 5: single stochastic trajectory
figure(5); clf; 
set(gcf, 'Color', 'w', 'Name', 'SSA Single Trajectory', 'NumberTitle', 'off');
hold on;

plot(T, Dyn_Conc(:,1), 'Color', [0, 0.5, 0], 'LineWidth', 0.5); 
yline(4.7, 'r--', 'LineWidth', 2, 'Label', 'Diseased State');
yline(1.45, 'b--', 'LineWidth', 2, 'Label', 'Healthy State');

xlabel('Time (t)');
ylabel('Amyloid Concentration (x)');
title(['Stochastic Simulation (SSA) - Volume \Omega = ' num2str(Omega)]);
legend({'Gillespie Trajectory', 'Diseased Attractor', 'Healthy Attractor'}, 'Location', 'Best');
grid on;
ylim([0 7]); 

fprintf('Final Amyloid Concentration: %.2f\n', Dyn_Conc(end,1));
if Dyn_Conc(end,1) > 3
    fprintf('RESULT: SWITCH OCCURRED.\n');
else
    fprintf('RESULT: NO SWITCH.\n');
end

end


%% --------------------------------------
%% PART 2: Intrinsic Amyloid Noise Analysis
%% --------------------------------------
function AD_DualNoise_Reproduction()
% The code keeps the Calcium volume fixed and loops through a list of decreasing Amyloid volumes

% Fixed parameters
V1 = 0.25;
V2 = 0.11;
V3 = 3;        
k1 = 0.3;
k2 = 5.0;
k3 = 1.0;
k4 = 1.0;
n  = 2;
epsilon = 0.01;

Omega_y = 25;
Omega_x_list = [25 20 15];

tMax = 4000;
burnFrac = 0.5;   
Nreal = 50;       
Conc_Start = [1.45, 0.35];

% Stochiometry
vMinus = [0 0; 0 0; 1 0; 0 0; 0 0; 0 1];
vPlus = [1 0; 1 0; 0 0; 0 1; 0 1; 0 0];

% FIGURE 6: Histograms
figure(6); clf; 
set(gcf,'Color','w', 'Name', 'SSA Dual Noise Histograms', 'NumberTitle', 'off');

% Loop over Omega_x
for idx = 1:length(Omega_x_list)
    
    Omega_x = Omega_x_list(idx);
    x_all = [];   
    
    fprintf('\nOmega_x = %d\n',Omega_x)
    
    c = [
        epsilon * V1 * Omega_x;
        0;
        epsilon * k1;
        V2 * Omega_y;
        k4;
        k2
    ];
    
    HillParams.V3 = V3;
    HillParams.K3 = k3;
    HillParams.n  = n;
    HillParams.epsilon = epsilon;
    HillParams.Omega_x = Omega_x;
    HillParams.Omega_y = Omega_y;
    
    parfor r = 1:Nreal
        
        InitialState = [
            round(Conc_Start(1) * Omega_x), ...
            round(Conc_Start(2) * Omega_y)
        ];
        
        [~, Dyn] = simDM_AD_Variant1(vMinus, vPlus, c, InitialState, tMax, HillParams);
        
        Dyn_Conc = [
            Dyn(:,1) / Omega_x, ...
            Dyn(:,2) / Omega_y
        ];
        
        burnin = floor(burnFrac * length(Dyn_Conc));
        x_ss = Dyn_Conc(burnin:end,1);
        
        x_all = [x_all; x_ss];
    end
    
    subplot(1,3,idx)
    histogram(x_all, 60, 'Normalization','pdf', ...
        'FaceColor',[0.2 0.8 0.2], 'EdgeColor','none')
    
    xlabel('x')
    ylabel('p_s(x)')
    title(['\Omega_x = ' num2str(Omega_x)])
    grid on
    xlim([0 10])
    
    fprintf('mean = %.3f, std = %.3f\n', mean(x_all), std(x_all))
end

sgtitle('Ensemble-averaged stationary distributions p_s(x)')
end

%% --------------------------------------
%% PART 3: Instrinsic Calcium Noise Analysis
%% --------------------------------------
function AD_Reproduction_Fig16()
% P-bifurcation: mean amyloid level vs calcium noise intensity.

% Fixed parameters
V1 = 0.25;
V2 = 0.11;
V3 = 3.0;
k1 = 0.30;
k2 = 5.0;
k3 = 1.0;
k4 = 1.0;
n  = 2;
epsilon = 0.01;

Omega_x = 20;                     
Omega_y_list = [50 20 15 10 5 2];  

tMax = 5000;
burnFrac = 0.5;
Nreal = 40;     
Conc_Start = [1.45, 0.35];

vMinus = [0 0; 0 0; 1 0; 0 0; 0 0; 0 1];
vPlus  = [1 0; 1 0; 0 0; 0 1; 0 1; 0 0];

mean_x = zeros(length(Omega_y_list), 1);

% Loop over Omega_y
for j = 1:length(Omega_y_list)
    
    Omega_y = Omega_y_list(j);
    noise_val = 1/sqrt(Omega_y);
    x_means = zeros(Nreal, 1);
    
    fprintf('Simulating Omega_y = %d (Noise Intensity = %.3f)...\n', Omega_y, noise_val)
    
    c = [
        epsilon * V1 * Omega_x;
        0;                      
        epsilon * k1;
        V2 * Omega_y;
        k4;
        k2
    ];
    
    HillParams.V3 = V3;
    HillParams.K3 = k3;
    HillParams.n  = n;
    HillParams.epsilon = epsilon;
    HillParams.Omega_x = Omega_x;
    HillParams.Omega_y = Omega_y;
    
    parfor r = 1:Nreal   
        InitialState = [
            round(Conc_Start(1) * Omega_x), ...
            round(Conc_Start(2) * Omega_y)
        ];
        
        [~, Dyn] = simDM_AD_Variant1(vMinus, vPlus, c, InitialState, tMax, HillParams);
        
        Dyn_Conc = [
            Dyn(:,1) / Omega_x, ...
            Dyn(:,2) / Omega_y
        ];
        
        burnin = floor(burnFrac * length(Dyn_Conc));
        if burnin >= length(Dyn_Conc), burnin = 1; end 
        x_means(r) = mean(Dyn_Conc(burnin:end, 1));
    end
    
    mean_x(j) = mean(x_means);
end

noise_intensity = 1 ./ sqrt(Omega_y_list);

% FIGURE 7: P-Bifurcation
figure(7); clf;
set(gcf, 'Color', 'w', 'Name', 'SSA P-Bifurcation', 'NumberTitle', 'off');

plot(noise_intensity, mean_x, 'bo-', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'MarkerSize', 8)

xlabel('Calcium Noise Intensity (\approx 1/\surd\Omega_y)')
ylabel('Average Amyloid Concentration')
title({'Mean Amyloid vs Calcium Noise'})
grid on

end


%% --------------------------------------
%% HELPER FUNCTIONS (DIRECT METHOD IMPLEMENTATION)
%% --------------------------------------

% VARIANT 1: Extracts Omega_x and Omega_y from HillParams struct
function [T, Dynamics] = simDM_AD_Variant1(vMinus, vPlus, c, initialState, tMax, HillParams)
    max_steps = 2000000; 
    T = nan(max_steps, 1);
    Dynamics = nan(max_steps, 2);
    
    v = vPlus - vMinus;
    i = 1;
    T(i) = 0;
    Dynamics(i,:) = initialState;
    
    V3 = HillParams.V3;
    K3 = HillParams.K3;
    n  = HillParams.n;
    epsilon = HillParams.epsilon;
    Omega_x = HillParams.Omega_x;
    Omega_y = HillParams.Omega_y;
    
    while T(i) < tMax
        X = Dynamics(i,1);
        Y = Dynamics(i,2);
        
        % Concentration for Hill function
        Y_conc = Y / Omega_y;
        
        % Update propensities
        a = zeros(1,6);
        a(1) = c(1);
        a(2) = epsilon * V3 * Omega_x * (Y_conc^n)/(K3^n + Y_conc^n);
        a(3) = c(3) * X;
        a(4) = c(4);
        a(5) = c(5) * X;
        a(6) = c(6) * Y;
        
        a0 = sum(a);
        if a0 == 0 || i >= max_steps, break; end
        
        r = rand(1,2);
        
        % Gillespie Direct Method
        tau = (1/a0)*log(1/r(1));
        
        % Find reaction index mu
        target = r(2) * a0;
        current_sum = 0;
        mu = 1;
        for k = 1:6
            current_sum = current_sum + a(k);
            if current_sum >= target
                mu = k;
                break;
            end
        end
        
        i = i + 1;
        T(i) = T(i-1) + tau;
        Dynamics(i,:) = Dynamics(i-1,:) + v(mu,:);
    end
    
    T = T(1:i);
    Dynamics = Dynamics(1:i,:);
end

% VARIANT 2: Accepts 'Omega' as an explicit argument
function [T, Dynamics] = simDM_AD_Variant2(vMinus, vPlus, c, initialState, tMax, Omega, HillParams)
    initialLenght = 100000; 
    T = nan(initialLenght,1); 
    Dynamics = nan(initialLenght,length(initialState)); 
    
    v = vPlus - vMinus;
    i = 1;
    T(i) = 0;
    Dynamics(i,:) = initialState;
    
    V3 = HillParams.V3;
    K3 = HillParams.K3;
    n  = HillParams.n;
    epsilon = HillParams.epsilon;
    
    while T(i) < tMax
        X_molecules = Dynamics(i,1); 
        Y_molecules = Dynamics(i,2); 
        Y_conc = Y_molecules / Omega;
        
        a = zeros(1, 6);
        a(1) = c(1); 
        hill_term = (Y_conc^n) / (K3^n + Y_conc^n);
        a(2) = epsilon * V3 * Omega * hill_term; 
        a(3) = c(3) * X_molecules; 
        a(4) = c(4); 
        a(5) = c(5) * X_molecules; 
        a(6) = c(6) * Y_molecules; 
        
        a0 = sum(a);
        if a0 == 0, break; end 
        
        r = rand(1,2);
        a_cumsum = cumsum(a);
        mu = find(a_cumsum >= r(1) * a0, 1);
        tau = (1/a0)*log(1/r(2));

        i = i+1;
        if i > length(T)
            T = [T; nan(initialLenght,1)];
            Dynamics = [Dynamics; nan(initialLenght,length(initialState))];
        end
        
        if (T(i-1) + tau <= tMax)
            T(i) = T(i-1) + tau;
            Dynamics(i,:) = Dynamics(i-1,:) + v(mu,:); 
        else
            T(i) = tMax;
            Dynamics(i,:) = Dynamics(i-1,:);
        end
    end
    
    T = T(1:i);
    Dynamics = Dynamics(1:i,:);
end
