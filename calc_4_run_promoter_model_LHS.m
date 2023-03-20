close all; clearvars; clc;

%% Select parent folder

parent_folder = 'D:\Google Drive\light_sweep_shared';

if ~isfolder(parent_folder)
    parent_folder = 'D:\GoogleDriveUW\light_sweep_shared';
end

% Get parameters for Msn2 localization function
light_sweep_folder = fullfile(parent_folder,'light_sweep_experiments');
model_folder = fullfile(parent_folder,'promoter_model');

% Set output folder
output_folder = fullfile(model_folder,'output');

% Import bounds for LHS sampling
opts = detectImportOptions(fullfile(model_folder,'model_params.xlsx'),'Sheet','promoter_params_LHS_bounds');
promoter_params_bounds = readtable(fullfile(model_folder,'model_params.xlsx'),opts);
promoter_params_labels = promoter_params_bounds.Properties.VariableNames;

% Import combined Msn2 measurements
load(fullfile(light_sweep_folder,'data_Msn2.mat'))
data_Msn2 = data_Msn2(data_Msn2.condition<=14,:);

%% Set simulation parameters

condition_list = unique(data_Msn2.condition);
t_measured = unique(data_Msn2.time);

initial_conditions = [1 0 0 0 0 0];
Kd_scale = 1;
fraction_active = 1;

number_simulation_rounds = 20;
guesses_per_round = 5E3;
number_guesses = number_simulation_rounds*guesses_per_round;

% LHS sampling
promoter_params_bounds.d2(1) = 0.0001; %%% Manual decrease d1%%%
promoter_params_bounds_lower = log(promoter_params_bounds{1,:});
promoter_params_bounds_upper = log(promoter_params_bounds{2,:});

promoter_params_LHS = lhsu(promoter_params_bounds_lower,promoter_params_bounds_upper,number_guesses);
promoter_params_LHS = exp(promoter_params_LHS);
promoter_params_LHS = array2table(promoter_params_LHS,'VariableNames',promoter_params_labels);

% save(fullfile(output_folder,'promoter_params_LHS'),'promoter_params_LHS');
return
%% Get Msn2 values and plot

% Organize Msn2 vs time
Msn2_measured_all = zeros(numel(t_measured),numel(condition_list));
for condition = 1:numel(condition_list)
    Msn2_measured_all(:,condition) = data_Msn2.mScarlet_localization(data_Msn2.condition==condition);
end

% Generate grouping indices for plotting Msn2
group = repmat(1:14,size(Msn2_measured_all,1),1); 
group = reshape(group,size(Msn2_measured_all,1)*size(Msn2_measured_all,2),1);
time = repmat(t_measured,size(Msn2_measured_all,2),1);

clear g; clc; close all
figure('units','normalized','outerposition',[0 0 0.7 0.5]);
g = gramm('x',data_Msn2.time,'y',data_Msn2.mScarlet_localization);
g.facet_wrap(data_Msn2.condition,'ncols',7);
g.geom_line();
g.set_line_options('base_size',1);
g.set_color_options('map',[55 55 55]/255);
g.set_names('x','time (min.)','y','Nuclear Msn2','column','');

g.draw;


%% Run model
close all; clc

% Initialize output variable
mCitrine_out = zeros(guesses_per_round,1,numel(t_measured),numel(condition_list));

% Loop through multiple simulation rounds (breaks simulation into multiple
% rounds (defined by number_simulation_rounds and number_guess_per_round)
% in case it crashes or is terminated early
for simulation_round = 1:number_simulation_rounds
    
    % Get parameters from latin hypercube
    idx_start = (simulation_round - 1)*guesses_per_round + 1;
    idx_end = simulation_round*guesses_per_round;
    promoter_params_LHS_round = promoter_params_LHS(idx_start:idx_end,:);
    
    % Run model with parameters from LHS and compare to measurements
    for condition_idx = 1:numel(condition_list)
        condition = condition_list(condition_idx);
        disp(condition)
        
        % Get Msn2 values
        Msn2_measured_temp = Msn2_measured_all(:,condition);

        % Loop through promotor parameter guesses and simulate model output
        parfor guess = 1:guesses_per_round
            
            % Get promoter parameters guess
            promoter_params_temp = promoter_params_LHS_round{guess,:};
            
            % Get modeled mCitrine
            [~,y] = ode45(@(t,y) promoter_ODEs(t,y,t_measured,Msn2_measured_temp,promoter_params_temp,Kd_scale,fraction_active),...
                t_measured,initial_conditions);
            
            % Get/store predicted mYFP measurements
            mCitrine_model_temp = y(:,6);            
            mCitrine_model(guess,1,:,condition) = mCitrine_model_temp;
            
        end
        
    end
    
    % Save model results for current round
    output_filename = strcat("mCitrine_model_round_",sprintf('%02d',simulation_round),".mat");
    save(fullfile(output_folder,output_filename),'mCitrine_model')
    
end

