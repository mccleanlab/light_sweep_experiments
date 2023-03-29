close all; clearvars; clc;

%% %%%%%%%%%%%% LOAD PRECALCULATED VARIABLES (IF APPLICABLE) %%%%%%%%%%%%%

% load('plot_figures_1_data.mat')
% return

%% %%%%%%%%%%%% IMPORT DATA AND CALCULATE PROMOTER THRESHOLDS %%%%%%%%%%%%%
%% Import measurements, simulations, etc

%%%%%%%%%%%%%%%%%%%%%%% Select input/output folders %%%%%%%%%%%%%%%%%%%%%%%

parent_folder  = pwd;
if contains(parent_folder,'\')
    parent_folder_idx   = strfind(parent_folder,'\');
else
    parent_folder_idx   = strfind(parent_folder,'/');
end
parent_folder = parent_folder(1:parent_folder_idx(end)-1);

%%%%%%%%%%%%%%%%%%%%%%%% Import data and parameters %%%%%%%%%%%%%%%%%%%%%%%

% Import LHS promoter parameters
promoter_model_folder = fullfile(parent_folder,'promoter_model');
model_solutions_folder = fullfile(promoter_model_folder,'output_1_Hill_1_Msn2_100k_d2');
load(fullfile(model_solutions_folder,'promoter_params_LHS.mat'));

% Load measurements
load(fullfile(parent_folder,'light_sweep_experiments','data.mat'))
load(fullfile(parent_folder,'light_sweep_experiments','mCitrine_stats.mat'))
load(fullfile(parent_folder,'light_sweep_experiments','data_stats.mat'))
load(fullfile(parent_folder,'light_sweep_experiments','data_Msn2.mat'))
load(fullfile(parent_folder,'light_sweep_experiments','data_Msn2_AUC.mat'))
load(fullfile(parent_folder,'light_vs_localization_experiments','data_light_dose.mat'))
load(fullfile(parent_folder,'promoter_analysis','reporter_promoters.mat'))
load(fullfile(parent_folder,'light_sweep_experiments','reporters_STRE_count.mat'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Process %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_stats = data_stats(data_stats.condition<=14,:);
data_stats = data_stats(ismember(data_stats.plasmid,{'pMM0845','pMM0846','pMM0847','pMM1079','pMM1070'}),:);
data_stats.mCitrine_cell = data_stats.mCitrine_cell - data_stats.mCitrine_cell_basal;
reporters_STRE_count.reporter = categorical(reporters_STRE_count.reporter);

% Import parameters for ideal Msn2 localization function
opts = detectImportOptions(fullfile(parent_folder,'plot_settings.xlsx'),'Sheet','Msn2_CT_params');
Msn2_params_list = readtable(fullfile(parent_folder,'plot_settings.xlsx'),opts);

% Get plot colors and Msn2_tex labels
opts = detectImportOptions(fullfile(parent_folder,'plot_settings.xlsx'),'Sheet','Msn2_tex');
opts.VariableTypes = {'categorical','categorical','double','double','double'};
plot_colors = readtable(fullfile(parent_folder,'plot_settings.xlsx'),opts);
plot_colors.RGB = [plot_colors.R, plot_colors.G, plot_colors.B];
plot_colors = plot_colors(:,{'Msn2','Msn2_tex','RGB'});
plot_colors.RGB(plot_colors.Msn2=='Msn2',:) = [70, 115, 190];

%%%%%%%%%%%%%%%%%%%% Get lists of reporters, Msn2, etc %%%%%%%%%%%%%%%%%%%%
Msn2_to_plot = {'Msn2(WT|4E|A)','Msn2(WT|4E|WT)','Msn2(WT|4E|T)'};
param_list = promoter_params_LHS.Properties.VariableNames;

strain_list = unique(data_stats.strain,'stable');
[strain_list, strain_order] = sort(string(strain_list));
strain_list = categorical(strain_list);

reporter_list = unique(data_stats.reporter,'stable');
reporter_list = reporter_list(strain_order);
reporter_list = categorical(reporter_list);

plasmid_list = unique(data_stats.plasmid);
plasmid_list = categorical(sort(string(plasmid_list)));
plasmid_list = plasmid_list(plasmid_list~='pMM1070');

%%%%%%%%%%%%%%%%%%%%%%%%%%% Import model output %%%%%%%%%%%%%%%%%%%%%%%%%%%
data_store = fileDatastore(fullfile(model_solutions_folder,'mCitrine_model_round_*.mat'),'ReadFcn',@load);
mCitrine_model = cell(size(data_store.Files,1),1);
idx = 1;

close all
initial_conditions = [1 0 0 0 0 0];
n_guesses_plot = 10;
fraction_active = 1;
reporters_to_plot = {'glpT','TKL2','ALD3','CTT1','DCS2','HXK1','RTN2','SIP18','SIP18 D6','SIP18 A4','DDR2','HSP12'};

% Loop through files in datastore and load
while hasdata(data_store)
    data_temp = read(data_store);
    
    field_name = fieldnames(data_temp);
    data_temp = data_temp.(field_name{:});
    
    mCitrine_model{idx,1} = data_temp;
    idx = idx + 1;
end

% Combine objects from datastore and sort by RSS
mCitrine_model = cat(1,mCitrine_model{:});

%% %%%%%%%%%%%%%%%%%% Fit measurements to ODE solutions %%%%%%%%%%%%%%%%%%%
clc

% Set parameters
n_guesses = size(mCitrine_model,1);
t_measured = unique(data_stats.time);
condition_list = unique(data_stats.condition);

% Organize Msn2 vs time
Msn2_measured_all = zeros(numel(t_measured),numel(condition_list));
for condition = 1:numel(condition_list)
    Msn2_measured_all(:,condition) = data_Msn2.mScarlet_localization(data_Msn2.condition==condition);
end


% Calculate norm of RSS for given strain/plasmid
promoter_fits = cell(numel(strain_list)*numel(plasmid_list),1);
idx = 1;
for strain_idx = 1:numel(strain_list)
    strain = strain_list(strain_idx);
    reporter = unique(data_stats.reporter(data_stats.strain==strain));
    
    for plasmid_idx = 1:numel(plasmid_list)
        
        % Get labels for strain/plasmid
        plasmid = plasmid_list(plasmid_idx);
        Msn2 = unique(data_stats.Msn2(data_stats.plasmid==plasmid));
        Msn2_tex = unique(data_stats.Msn2_tex(data_stats.plasmid==plasmid));
        
        % Organize measurements for strain/plasmid
        replicate_list = unique(data_stats.replicate(data_stats.strain==strain & data_stats.plasmid==plasmid));
        mCitrine_measured = zeros(1,numel(replicate_list),numel(t_measured),numel(condition_list));
        
        for replicate_idx = 1:numel(replicate_list)
            replicate = replicate_list(replicate_idx);
            for condition_idx = 1:numel(condition_list)
                condition = condition_list(condition_idx);
                
                mCitrine_measured(1,replicate_idx,:,condition_idx) = data_stats.mCitrine_cell(...
                    data_stats.strain==strain & data_stats.plasmid==plasmid ...
                    & data_stats.condition==condition & data_stats.replicate==replicate);
            end
        end
        
        
        % Calculate norm of RSS
        RSS = nansum((mCitrine_measured - mCitrine_model).^2,2:3);
        RSS = squeeze(RSS);
        RSS_norm = sqrt(sum(RSS.^2,2));
        
        RSS_sp = nansum((mCitrine_measured(:,:,:,1:9) - mCitrine_model(:,:,:,1:9)).^2,2:3);
        RSS_sp = squeeze(RSS_sp);
        RSS_norm_sp = sqrt(sum(RSS_sp.^2,2));
        
        % Sort parameters by norm of RSS for given strain/plasmid
        [RSS_norm, RSS_sort_idx] = sort(RSS_norm,'ascend');
        RSS_norm_sp = RSS_norm_sp(RSS_sort_idx,:);
        promoter_params_sort_RSS = promoter_params_LHS(RSS_sort_idx,:);
        
        % Save best parameters
        promoter_fits_temp = table();
        promoter_fits_temp.strain(1:n_guesses,1) = strain;
        promoter_fits_temp.reporter(1:n_guesses,1) = reporter;
        promoter_fits_temp.plasmid(1:n_guesses,1) = plasmid;
        promoter_fits_temp.Msn2(1:n_guesses,1) = Msn2;
        promoter_fits_temp.Msn2_tex(1:n_guesses,1) = Msn2_tex;
        promoter_fits_temp = [promoter_fits_temp,promoter_params_sort_RSS(1:n_guesses,:)];
        promoter_fits_temp.RSS_norm(1:n_guesses,1) = RSS_norm(1:n_guesses,:);
        promoter_fits_temp.RSS_norm_fc(1:n_guesses,1) = RSS_norm(1:n_guesses,:)/min(RSS_norm);
        promoter_fits_temp.rank(1:n_guesses,1) = (1:n_guesses)';
        promoter_fits{idx,1} = promoter_fits_temp;
        idx = idx + 1;
        
    end
end

promoter_fits = vertcat(promoter_fits{:});
promoter_fits.Kd = promoter_fits.K.^promoter_fits.n;

%% Calculate promoter responses to measured pulses Msn2 

%%%%%%% Model promoter responses to measured pulses Msn2 w/ K scaling %%%%%

K_scale_list = [0.5, 1, 2];
promoter_response_measured = cell(numel(reporters_to_plot)*numel(plasmid_list)*n_guesses_plot*numel(K_scale_list)*14,21);
idx = 1;
for ii = 1:numel(reporters_to_plot)
    reporter = categorical(reporters_to_plot(ii));
    strain = unique(data_stats.strain(data_stats.reporter==reporter));
    disp(strain)
    
    for plasmid_idx = 1:numel(plasmid_list)
        
        plasmid = plasmid_list(plasmid_idx);
        Msn2 = unique(data_stats.Msn2(data_stats.plasmid==plasmid));
        Msn2_tex = unique(data_stats.Msn2_tex(data_stats.plasmid==plasmid));
        
        % Get top params for strain/plasmid
        subset = promoter_fits.strain==strain & promoter_fits.plasmid==plasmid;
        promoter_params_strain = promoter_fits{subset,param_list};
        promoter_params_strain = promoter_params_strain(1:n_guesses_plot,:);
        
        for guess = 1:n_guesses_plot
            
            promoter_params_temp = promoter_params_strain(guess,:);
            
            for K_scale_idx = 1:numel(K_scale_list)
                
                K_scale = K_scale_list(K_scale_idx);
                
                % Initialize variables
                P_unbound_model_temp = zeros(numel(t_measured),numel(condition_list));
                P_bound_model_temp = zeros(numel(t_measured),numel(condition_list));
                P_active_model_temp = zeros(numel(t_measured),numel(condition_list));
                mRNA_model_temp = zeros(numel(t_measured),numel(condition_list));
                mCitrine_model_temp = zeros(numel(t_measured),numel(condition_list));
                
                parfor condition = 1:14
                    Msn2_measured_temp = Msn2_measured_all(:,condition);
                    
                    [~,y] = ode45(@(t,y) promoter_ODEs(t,y,t_measured,Msn2_measured_temp,promoter_params_temp,K_scale,fraction_active),...
                        t_measured,initial_conditions);
                    P_unbound_model_temp(:,condition) = y(:,1);
                    P_bound_model_temp(:,condition) = y(:,2);
                    P_active_model_temp(:,condition) = y(:,3);
                    mRNA_model_temp(:,condition) = y(:,4);
                    mCitrine_model_temp(:,condition) = y(:,6);
                end
                
                for condition = 1:14
                    condition_display = unique(data_Msn2.condition_display(data_Msn2.condition==condition));   

                    promoter_response_measured{idx,1} = strain;
                    promoter_response_measured{idx,2} = reporter;
                    promoter_response_measured{idx,3} = plasmid;
                    promoter_response_measured{idx,4} = Msn2;
                    promoter_response_measured{idx,5} = Msn2_tex;
                    
                    promoter_response_measured{idx,6} = condition;
                     promoter_response_measured{idx,7} = condition_display;
                    promoter_response_measured{idx,8} = t_measured;
                    promoter_response_measured{idx,9} = Msn2_measured_all(:,condition);
                    promoter_response_measured{idx,10} = trapz(t_measured,Msn2_measured_all(:,condition));
                    
                    promoter_response_measured{idx,11} = guess;
                    promoter_response_measured{idx,12} = promoter_params_temp;
                    promoter_response_measured{idx,13} = K_scale;
                    promoter_response_measured{idx,14} = fraction_active;
                    
                    promoter_response_measured{idx,15} = P_unbound_model_temp(:,condition);
                    promoter_response_measured{idx,16} = P_bound_model_temp(:,condition);
                    promoter_response_measured{idx,17} = P_active_model_temp(:,condition);
                    promoter_response_measured{idx,18} = mRNA_model_temp(:,condition);
                    promoter_response_measured{idx,19} = mCitrine_model_temp(:,condition);
                    promoter_response_measured{idx,20} = max(P_bound_model_temp(:,condition));
                    promoter_response_measured{idx,21} = max(mCitrine_model_temp(:,condition));
                    
                    idx = idx + 1;
                end
            end
        end
    end
end

promoter_response_measured = cell2table(promoter_response_measured,'VariableNames',...
    {'strain','reporter','plasmid','Msn2','Msn2_tex',...
    'condition','condition_display','time','mScarlet_nuclear','mScarlet_AUC',...
    'guess','params','K_scale','fraction_active',...
    'P_unbound','P_bound','P_active','mRNA','mCitrine','P_active_max','mCitrine_max'});

promoter_response_measured = splitvars(promoter_response_measured,'params','NewVariableNames',param_list);

%% Calculate promoter response to ideal pulses of nuclear Msn2 

%%%%%%%%%%%%%%%%%%%% Define ideal pulses nuclear Msn2 %%%%%%%%%%%%%%%%%%%%%
pulse_end_list = [0:2.5:7.5,10:5:75];
pulse_amplitude_list = 0:0.10:1.5;
pulse_t = linspace(min(t_measured),max(t_measured),1000)';

Msn2_pulses_ideal = cell(numel(pulse_end_list) + numel(pulse_amplitude_list),11);
idx = 1;
for ii = 1:numel(pulse_end_list)
    
    t1 = pulse_end_list(ii);
    A = 1;
    
    Msn2_params = Msn2_params_list(Msn2_params_list.condition==9,:);
    Msn2_params.t1 = pulse_end_list(ii);
    Msn2_params.A = 1;
    
    Msn2_pulses_ideal{idx,1} = idx;
    Msn2_pulses_ideal{idx,2} = Msn2_params.signal_type;
    Msn2_pulses_ideal{idx,3} = Msn2_params.A;
    Msn2_pulses_ideal{idx,4} = Msn2_params.t0;
    Msn2_pulses_ideal{idx,5} = Msn2_params.t1;
    Msn2_pulses_ideal{idx,6} = Msn2_params.t2;
    Msn2_pulses_ideal{idx,7} = Msn2_params.cycles;
    Msn2_pulses_ideal{idx,8} = Msn2_params.c1;
    Msn2_pulses_ideal{idx,9} = Msn2_params.c2;
    Msn2_pulses_ideal{idx,10} = categorical("duration");
    Msn2_pulses_ideal{idx,11} = pulse_t;
    Msn2_pulses_ideal{idx,12} = Msn2_CT(pulse_t,Msn2_params);
    Msn2_pulses_ideal{idx,13} = trapz(pulse_t,Msn2_CT(pulse_t,Msn2_params));
    
    idx = idx + 1;
end

for jj = 1:numel(pulse_amplitude_list)
    
    Msn2_params = Msn2_params_list(Msn2_params_list.condition==9,:);
    Msn2_params.t1 = 50;
    Msn2_params.A = pulse_amplitude_list(jj);
    
    Msn2_pulses_ideal{idx,1} = idx;
    Msn2_pulses_ideal{idx,2} = Msn2_params.signal_type;
    Msn2_pulses_ideal{idx,3} = Msn2_params.A;
    Msn2_pulses_ideal{idx,4} = Msn2_params.t0;
    Msn2_pulses_ideal{idx,5} = Msn2_params.t1;
    Msn2_pulses_ideal{idx,6} = Msn2_params.t2;
    Msn2_pulses_ideal{idx,7} = Msn2_params.cycles;
    Msn2_pulses_ideal{idx,8} = Msn2_params.c1;
    Msn2_pulses_ideal{idx,9} = Msn2_params.c2;
    Msn2_pulses_ideal{idx,10} = categorical("amplitude");
    Msn2_pulses_ideal{idx,11} = pulse_t;
    Msn2_pulses_ideal{idx,12} = Msn2_CT(pulse_t,Msn2_params);
    Msn2_pulses_ideal{idx,13} = trapz(pulse_t,Msn2_CT(pulse_t,Msn2_params));
    
    idx = idx + 1;
end

for condition = 10:14
    
    Msn2_params = Msn2_params_list(Msn2_params_list.condition==condition,:);
    pulse_y = Msn2_CT(pulse_t,Msn2_params);
    
    Msn2_pulses_ideal{idx,1} = idx;
    Msn2_pulses_ideal{idx,2} = Msn2_params.signal_type;
    Msn2_pulses_ideal{idx,3} = Msn2_params.A;
    Msn2_pulses_ideal{idx,4} = Msn2_params.t0;
    Msn2_pulses_ideal{idx,5} = Msn2_params.t1;
    Msn2_pulses_ideal{idx,6} = Msn2_params.t2;
    Msn2_pulses_ideal{idx,7} = Msn2_params.cycles;
    Msn2_pulses_ideal{idx,8} = Msn2_params.c1;
    Msn2_pulses_ideal{idx,9} = Msn2_params.c2;
    Msn2_pulses_ideal{idx,10} = categorical("pulsed");
    Msn2_pulses_ideal{idx,11} = pulse_t;
    Msn2_pulses_ideal{idx,12} = Msn2_CT(pulse_t,Msn2_params);
    Msn2_pulses_ideal{idx,13} = trapz(pulse_t,Msn2_CT(pulse_t,Msn2_params));
    
    idx = idx + 1;
end

Msn2_pulses_ideal = cell2table(Msn2_pulses_ideal,'VariableNames',...
    {'pulse_idx','signal_type','A','t0','t1','t2','cycles','c1','c2','pulse_label','pulse_t','pulse_y','mScarlet_AUC'});
Msn2_pulses_ideal.pulse_idx(:,1) = 1:size(Msn2_pulses_ideal,1);
pulse_idx_list = unique(Msn2_pulses_ideal.pulse_idx);


%%%%%%%%%%%%%%%%%%% Calculate AUC of ideal Msn2 traces %%%%%%%%%%%%%%%%%%%%
t_Msn2_ideal = linspace(min(data_Msn2.time),max(data_Msn2.time),1000)';
t_measured = unique(data_Msn2.time);

idx = 1;
Msn2_ideal = cell(14,6);
for condition = 1:14
    Msn2_params = Msn2_params_list(Msn2_params_list.condition==condition,:);
    Msn2_localization_ideal = Msn2_CT(t_Msn2_ideal,Msn2_params);
    Msn2_AUC_ideal = trapz(t_Msn2_ideal,Msn2_localization_ideal);
    
    y_measured = data_Msn2.mScarlet_localization(data_Msn2.condition==condition);
    y_ideal = interp1(t_Msn2_ideal,Msn2_localization_ideal,t_measured);
    y_diff = y_ideal - y_measured;
    
    % Manually exlcude peaks where predicted and ideal offset by 1 frame
    exclude = y_diff<-0.025 | (t_measured<14 & y_diff>0.15);
    y_diff_exclude = y_diff;
    y_diff_exclude(exclude) = nan;
    
    Msn2_ideal{idx,1} = condition;
    Msn2_ideal{idx,2} = t_measured;
    Msn2_ideal{idx,3} = y_ideal;
    Msn2_ideal{idx,4} = y_diff;
    Msn2_ideal{idx,5} = y_diff_exclude;
    Msn2_ideal{idx,6} = Msn2_AUC_ideal;

    idx = idx + 1;
        
end

Msn2_ideal = cell2table(Msn2_ideal,'VariableNames',{'condition','time','Msn2_localization_ideal',...
    'Msn2_localization_diff','Msn2_localization_exclude','Msn2_AUC_ideal'});

% Join ideal traces back into measurements table
data_Msn2_AUC = outerjoin(data_Msn2_AUC,Msn2_ideal, 'Type', 'Left', 'MergeKeys', true);

%%%%%%%% Model promoter responses to ideal pulses Msn2 w/ K scaling %%%%%%%

promoter_response_ideal = cell(numel(reporter_list)*numel(plasmid_list)*n_guesses_plot*numel(K_scale_list)*size(Msn2_pulses_ideal,1),26);

idx = 1;
for reporter_idx = 1:numel(reporter_list)
    reporter = categorical(reporter_list(reporter_idx));
    strain = unique(data_stats.strain(data_stats.reporter==reporter));
    disp(strain)
    
    % Plot fits over measurements
    for plasmid_idx = 1:numel(plasmid_list)
        
        plasmid = plasmid_list(plasmid_idx);
        Msn2 = unique(data_stats.Msn2(data_stats.plasmid==plasmid));
        Msn2_tex = unique(data_stats.Msn2_tex(data_stats.plasmid==plasmid));
        
        % Get params for strain/plasmid
        subset = promoter_fits.strain==strain & promoter_fits.plasmid==plasmid;
        promoter_params_strain = promoter_fits(subset,param_list);
        promoter_params_strain = promoter_params_strain{1:n_guesses_plot,:};
        
        for guess = 1:n_guesses_plot
            
            promoter_params_temp = promoter_params_strain(guess,:);
            
            for K_scale_idx = 1:numel(K_scale_list)
                
                K_scale = K_scale_list(K_scale_idx);
                
                % Initialize variables
                P_unbound_model_temp = zeros(numel(t_measured),size(Msn2_pulses_ideal,1));
                P_bound_model_temp = zeros(numel(t_measured),size(Msn2_pulses_ideal,1));
                P_active_model_temp = zeros(numel(t_measured),size(Msn2_pulses_ideal,1));
                mRNA_model_temp = zeros(numel(t_measured),size(Msn2_pulses_ideal,1));
                mCitrine_model_temp = zeros(numel(t_measured),size(Msn2_pulses_ideal,1));
                
                parfor pulse_idx = 1:numel(pulse_idx_list)
                    t_temp = Msn2_pulses_ideal.pulse_t{Msn2_pulses_ideal.pulse_idx==pulse_idx};
                    Msn2_temp = Msn2_pulses_ideal.pulse_y{Msn2_pulses_ideal.pulse_idx==pulse_idx};
                    
                    [~,y] = ode45(@(t,y) promoter_ODEs(t,y,t_temp,Msn2_temp,promoter_params_temp,K_scale,fraction_active),...
                        t_measured,initial_conditions);
                    P_unbound_model_temp(:,pulse_idx) = y(:,1);
                    P_bound_model_temp(:,pulse_idx) = y(:,2);
                    P_active_model_temp(:,pulse_idx) = y(:,3);
                    mRNA_model_temp(:,pulse_idx) = y(:,4);
                    mCitrine_model_temp(:,pulse_idx) = y(:,6);
                end
                
                for pulse_idx = 1:numel(pulse_idx_list)
                    
                    t_temp = Msn2_pulses_ideal.pulse_t{Msn2_pulses_ideal.pulse_idx==pulse_idx};
                    Msn2_temp = Msn2_pulses_ideal.pulse_y{Msn2_pulses_ideal.pulse_idx==pulse_idx};
                    
                    pulse_label = Msn2_pulses_ideal.pulse_label(Msn2_pulses_ideal.pulse_idx==pulse_idx);
                    A = Msn2_pulses_ideal.A(Msn2_pulses_ideal.pulse_idx==pulse_idx);
                    t0 = Msn2_pulses_ideal.t0(Msn2_pulses_ideal.pulse_idx==pulse_idx);
                    t1 = Msn2_pulses_ideal.t1(Msn2_pulses_ideal.pulse_idx==pulse_idx);
                    t2 = Msn2_pulses_ideal.t2(Msn2_pulses_ideal.pulse_idx==pulse_idx);
                    cycles = Msn2_pulses_ideal.cycles(Msn2_pulses_ideal.pulse_idx==pulse_idx);
                    
                    promoter_response_ideal{idx,1} = strain;
                    promoter_response_ideal{idx,2} = reporter;
                    promoter_response_ideal{idx,3} = plasmid;
                    promoter_response_ideal{idx,4} = Msn2;
                    promoter_response_ideal{idx,5} = Msn2_tex;
                    
                    promoter_response_ideal{idx,6} = pulse_idx;
                    promoter_response_ideal{idx,7} = pulse_label;
                    promoter_response_ideal{idx,8} = A;
                    promoter_response_ideal{idx,9} = t0;
                    promoter_response_ideal{idx,10} = t1;
                    promoter_response_ideal{idx,11} = t2;
                    promoter_response_ideal{idx,12} = cycles;
                    
                    promoter_response_ideal{idx,13} = t_temp;
                    promoter_response_ideal{idx,14} = Msn2_temp;
                    promoter_response_ideal{idx,15} = trapz(t_temp,Msn2_temp);
                    
                    promoter_response_ideal{idx,13} = t_measured;
                    promoter_response_ideal{idx,14} = interp1(t_temp,Msn2_temp,t_measured);
                    promoter_response_ideal{idx,15} = trapz(t_temp,Msn2_temp);
                    
                    promoter_response_ideal{idx,16} = guess;
                    promoter_response_ideal{idx,17} = promoter_params_temp;
                    promoter_response_ideal{idx,18} = K_scale;
                    promoter_response_ideal{idx,19} = fraction_active;
                    
                    promoter_response_ideal{idx,20} = P_unbound_model_temp(:,pulse_idx);
                    promoter_response_ideal{idx,21} = P_bound_model_temp(:,pulse_idx);
                    promoter_response_ideal{idx,22} = P_active_model_temp(:,pulse_idx);
                    promoter_response_ideal{idx,23} = mRNA_model_temp(:,pulse_idx);
                    promoter_response_ideal{idx,24} = mCitrine_model_temp(:,pulse_idx);
                    promoter_response_ideal{idx,25} = max(P_active_model_temp(:,pulse_idx));
                    promoter_response_ideal{idx,26} = max(mCitrine_model_temp(:,pulse_idx));
                    
                    idx = idx + 1;
                    
                end
            end
        end
    end
end

promoter_response_ideal = cell2table(promoter_response_ideal,'VariableNames',...
    {'strain','reporter','plasmid','Msn2','Msn2_tex',...
    'pulse_idx','pulse_label','A','t0' 't1','t2','cycles',...
    'time','mScarlet_localization','mScarlet_AUC',...
    'guess','params','K_scale','fraction_active',...
    'P_unbound','P_bound','P_active','mRNA','mCitrine','P_active_max','mCitrine_max'});

% Calculate k3*P_active & k3*P_active_max
promoter_response_ideal = splitvars(promoter_response_ideal,'params','NewVariableNames',param_list);
promoter_response_ideal.P_active_max_k3 = promoter_response_ideal.P_active_max.*promoter_response_ideal.k3;
for ii = 1:size(promoter_response_ideal,1)
    promoter_response_ideal.P_active_k3(ii) = {promoter_response_ideal.P_active{ii}.*promoter_response_ideal.k3(ii)};
end


%% Calculate promoter thresholds from measured nuclear Msn2

% Initialize
promoter_thresholds_measured = cell(numel(strain_list)*numel(plasmid_list),1);

%%%%%%%%%%%%%% Calculate promoter thresholds per Msn2 mutant %%%%%%%%%%%%%%
idx = 1;
for ii = 1:numel(strain_list)
    strain = strain_list(ii);
    reporter = unique(mCitrine_stats.reporter(mCitrine_stats.strain==strain));
    
    for jj = 1:numel(plasmid_list)
        
        plasmid = plasmid_list(jj);
        Msn2 = unique(mCitrine_stats.Msn2(mCitrine_stats.plasmid==plasmid));
        Msn2_tex = unique(mCitrine_stats.Msn2_tex(mCitrine_stats.plasmid==plasmid));
        grp_vars = {'strain','plasmid','condition','amplitudes'};
        
        % Calculate amplitude threshold 
        conditions_to_fit = 1:2:9;
        subset = (mCitrine_stats.strain==strain & mCitrine_stats.plasmid==plasmid ...
            & ismember(mCitrine_stats.condition,conditions_to_fit));
        mCitrine_temp = mCitrine_stats(subset,:);
        
        replicate_list = unique(mCitrine_temp.replicate);
        amplitude_threshold = zeros(numel(replicate_list),1);
        amplitude_threshold_y = zeros(numel(replicate_list),1);
        for kk = 1:numel(replicate_list)
            replicate = replicate_list(kk);
            
            mCitrine_temp_replicate = mCitrine_temp(mCitrine_temp.replicate==replicate,:);
            x = linspace(min(mCitrine_temp_replicate.mScarlet_AUC),max(mCitrine_temp_replicate.mScarlet_AUC),1000);
            y = interp1(mCitrine_temp_replicate.mScarlet_AUC,mCitrine_temp_replicate.mCitrine_max,x);
            
            y_50 = 0.5*(max(y) - min(y)) + min(y);
            
            [~,x_50] = min(abs(y - y_50));
            amplitude_threshold(kk) = x(x_50);
            amplitude_threshold_y(kk) = y_50;
            
        end
        
        amplitude_threshold_std = std(amplitude_threshold);
        amplitude_threshold = mean(amplitude_threshold);
        amplitude_threshold_y = mean(amplitude_threshold_y);
        
        % Calculate duration threshold 
        conditions_to_fit = [1,2,4,6,8,9];
        subset = (mCitrine_stats.strain==strain & mCitrine_stats.plasmid==plasmid ...
            & ismember(mCitrine_stats.condition,conditions_to_fit));
        mCitrine_temp = mCitrine_stats(subset,:);
        
        replicate_list = unique(mCitrine_temp.replicate);
        duration_threshold = zeros(numel(replicate_list),1);
        duration_threshold_y = zeros(numel(replicate_list),1);
        for kk = 1:numel(replicate_list)
            
            replicate = replicate_list(kk);
            mCitrine_temp_replicate = mCitrine_temp(mCitrine_temp.replicate==replicate,:);
            x = linspace(min(mCitrine_temp_replicate.mScarlet_AUC),max(mCitrine_temp_replicate.mScarlet_AUC),1000);
            y = interp1(mCitrine_temp_replicate.mScarlet_AUC,mCitrine_temp_replicate.mCitrine_max,x);
            
            y_50 = 0.5*(max(y) - min(y)) + min(y);
            [~,x_50] = min(abs(y - y_50));
            duration_threshold(kk) = x(x_50);
            duration_threshold_y(kk) = y_50;
        end
        
        duration_threshold_std = std(duration_threshold);
        duration_threshold = mean(duration_threshold);
        duration_threshold_y = mean(duration_threshold_y);
        
       % Calculate slope ratio
        
        % Calculate slope of continuous conditions
        conditions_to_fit = [1,2,6,9];
        
        subset = (mCitrine_stats.strain==strain & mCitrine_stats.plasmid==plasmid ...
            & ismember(mCitrine_stats.condition,conditions_to_fit));
        mCitrine_temp = mCitrine_stats(subset,:);
        
        mdl_continuous = fitglm(mCitrine_temp.mScarlet_AUC,mCitrine_temp.mCitrine_max);
        m_continuous = mdl_continuous.Coefficients.Estimate(2);
        m_continuous_se = mdl_continuous.Coefficients.SE(2);
        m_continuous_p = mdl_continuous.Coefficients.pValue(2);
        b_continuous = mdl_continuous.Coefficients.Estimate(1);
        
        x_mdl_continuous = (linspace(min(mCitrine_temp.mScarlet_AUC),max(mCitrine_temp.mScarlet_AUC),1000))';
        y_mdl_continuous = m_continuous.*x_mdl_continuous + b_continuous;
        
        % Calculate slope of pulsed conditions
        conditions_to_fit = [1,10,11,12];
        
        subset = (mCitrine_stats.strain==strain & mCitrine_stats.plasmid==plasmid ...
            & ismember(mCitrine_stats.condition,conditions_to_fit));
        mCitrine_temp = mCitrine_stats(subset,:);
        
        mdl_pulsed = fitglm(mCitrine_temp.mScarlet_AUC,mCitrine_temp.mCitrine_max);
        m_pulsed = mdl_pulsed.Coefficients.Estimate(2);
        m_pulsed_se = mdl_pulsed.Coefficients.SE(2);
        m_pulsed_p = mdl_pulsed.Coefficients.pValue(2);
        b_pulsed = mdl_pulsed.Coefficients.Estimate(1);
        
        x_mdl_pulsed = (linspace(min(mCitrine_temp.mScarlet_AUC),max(mCitrine_temp.mScarlet_AUC),1000))';
        y_mdl_pulsed = m_pulsed.*x_mdl_pulsed + b_pulsed;
        
        slope_ratio = m_pulsed./m_continuous;
        slope_ratio_se = (m_pulsed/m_continuous)* sqrt((m_pulsed_se/m_pulsed).^2 + (m_continuous_se/m_continuous).^2);
        
        %  Calculate interpulse decay
        conditions_to_fit = [6,11,13,14];
        subset = (mCitrine_stats.strain==strain & mCitrine_stats.plasmid==plasmid ...
            & ismember(mCitrine_stats.condition,conditions_to_fit));
        mCitrine_temp = mCitrine_stats(subset,:);
        
        %         [mdl_interpulse , gof_mdl_interpulse] = fit(mCitrine_temp.mScarlet_AUC,mCitrine_temp.mCitrine_max,'exp1','StartPoint',[1,0.1]);
        interpulse_function = @(a,x) a(1)*exp(a(2)*x);
        %         interpulse_function = @(a,x) exp(a(1)*x) + a(2);
        mdl_interpulse = fitnlm(mCitrine_temp.mScarlet_AUC,mCitrine_temp.mCitrine_max,interpulse_function,[1,0.1]);
        
        b_interpulse = mdl_interpulse.Coefficients.Estimate(2);
        b_interpulse_se = mdl_interpulse.Coefficients.SE(2);
        b_interpulse_p = mdl_interpulse.Coefficients.pValue(2);
        a_interpulse = mdl_interpulse.Coefficients.Estimate(1);
        
        x_mdl_interpulse = (linspace(min(mCitrine_temp.mScarlet_AUC),max(mCitrine_temp.mScarlet_AUC),1000))';
        y_mdl_interpulse = a_interpulse*exp(b_interpulse.*x_mdl_interpulse);
        %         y_mdl_interpulse = exp(a_interpulse.*x_mdl_interpulse) + b_interpulse;
                
        % Collect measurements and thresholds
        promoter_thresholds_measured{idx,1} = strain;
        promoter_thresholds_measured{idx,2} = reporter;
        promoter_thresholds_measured{idx,3} = plasmid;
        promoter_thresholds_measured{idx,4} = Msn2;
        promoter_thresholds_measured{idx,5} = Msn2_tex;
        
        promoter_thresholds_measured{idx,6} = amplitude_threshold;
        promoter_thresholds_measured{idx,7} = amplitude_threshold_std;
        promoter_thresholds_measured{idx,8} = amplitude_threshold_y;
        
        promoter_thresholds_measured{idx,9} = duration_threshold;
        promoter_thresholds_measured{idx,10} = duration_threshold_std;
        promoter_thresholds_measured{idx,11} = duration_threshold_y;
        
        promoter_thresholds_measured{idx,12} = slope_ratio;
        promoter_thresholds_measured{idx,13} = slope_ratio_se;
        
        promoter_thresholds_measured{idx,14} = mdl_continuous;
        promoter_thresholds_measured{idx,15} = x_mdl_continuous;
        promoter_thresholds_measured{idx,16} = y_mdl_continuous;
        promoter_thresholds_measured{idx,17} = m_continuous;
        promoter_thresholds_measured{idx,18} = m_pulsed_se;
        promoter_thresholds_measured{idx,19} = m_continuous_p;
        
        promoter_thresholds_measured{idx,20} = mdl_pulsed;
        promoter_thresholds_measured{idx,21} = x_mdl_pulsed;
        promoter_thresholds_measured{idx,22} = y_mdl_pulsed;
        promoter_thresholds_measured{idx,23} = m_pulsed;
        promoter_thresholds_measured{idx,24} = m_pulsed_se;
        promoter_thresholds_measured{idx,25} = m_pulsed_p;
        
        promoter_thresholds_measured{idx,26} = mdl_interpulse;
        promoter_thresholds_measured{idx,27} = x_mdl_interpulse;
        promoter_thresholds_measured{idx,28} = y_mdl_interpulse;
        promoter_thresholds_measured{idx,29} = b_interpulse;
        promoter_thresholds_measured{idx,30} = b_interpulse_se;
        promoter_thresholds_measured{idx,31} = b_interpulse_p;
        
        idx = idx +1;
    end
end

promoter_thresholds_measured = cell2table(promoter_thresholds_measured,'VariableNames',...
    {'strain','reporter','plasmid','Msn2','Msn2_tex',...
    'amplitude_threshold','amplitude_threshold_std','amplitude_threshold_y',...
    'duration_threshold','duration_threshold_std','duration_threshold_y',...
    'slope_ratio','slope_ratio_se',...
    'mdl_slope_continuous','x_continuous','y_continuous','slope_continuous','slope_continuous_se','slope_continuous_p',...
    'mdl_slope_pulsed','x_pulsed','y_pulsed','slope_pulsed','slope_pulsed_se','slope_pulsed_p',...
    'mdl_slope_interpulse','x_interpulse','y_interpulse','interpulse_decay','interpulse_decay_se','interpulse_decay_p',...
    });

%%%%%%%%%%%%%%%%%%%%%% Calculate thresholds vs Msn2* %%%%%%%%%%%%%%%%%%%%%%
promoter_thresholds_measured.amplitude_threshold_vs_WT(:,1) = nan;
promoter_thresholds_measured.amplitude_threshold_vs_WT_std(:,1) = nan;
promoter_thresholds_measured.duration_threshold_vs_WT(:,1) = nan;
promoter_thresholds_measured.duration_threshold_vs_WT_std(:,1) = nan;

for ii = 1:numel(strain_list)
    
    strain = strain_list(ii);
    reporter = unique(mCitrine_stats.reporter(mCitrine_stats.strain==strain));
    
    amplitude_threshold_y_pMM0847 = promoter_thresholds_measured.amplitude_threshold_y(...
        promoter_thresholds_measured.strain==strain & promoter_thresholds_measured.plasmid=='pMM0847',:);
    
    duration_threshold_y_pMM0847 = promoter_thresholds_measured.duration_threshold_y(...
        promoter_thresholds_measured.strain==strain & promoter_thresholds_measured.plasmid=='pMM0847',:);
    
    for jj = 1:numel(plasmid_list)
        plasmid = plasmid_list(jj);
        
        % Calculate y value for WT amplitude threshold
        conditions_to_fit = 1:2:9;
        subset = (mCitrine_stats.strain==strain & mCitrine_stats.plasmid==plasmid ...
            & ismember(mCitrine_stats.condition,conditions_to_fit));
        mCitrine_temp = mCitrine_stats(subset,:);
        
        replicate_list = unique(mCitrine_temp.replicate);
        amplitude_threshold_vs_WT = zeros(numel(replicate_list),1);
        
        % Calculate amplitude threshold vs WT
        y_diff_amplitude = zeros(numel(replicate_list),1);
        for kk = 1:numel(replicate_list)
            replicate = replicate_list(kk);
            
            mCitrine_temp_replicate = mCitrine_temp(mCitrine_temp.replicate==replicate,:);
            x = linspace(min(mCitrine_temp_replicate.mScarlet_AUC),max(mCitrine_temp_replicate.mScarlet_AUC),1000);
            y = interp1(mCitrine_temp_replicate.mScarlet_AUC,mCitrine_temp_replicate.mCitrine_max,x);
            
            [y_diff,x_50] = min(abs(y - amplitude_threshold_y_pMM0847));
            y_diff_amplitude(kk) = y_diff;
            amplitude_threshold_vs_WT(kk) = x(x_50);
            
        end
        
        amplitude_threshold_vs_WT_std = std(amplitude_threshold_vs_WT);
        amplitude_threshold_vs_WT = mean(amplitude_threshold_vs_WT);
        
        % Calculate y value for WT duration threshold
        conditions_to_fit = [1,2,4,6,8,9];
        subset = (mCitrine_stats.strain==strain & mCitrine_stats.plasmid==plasmid ...
            & ismember(mCitrine_stats.condition,conditions_to_fit));
        mCitrine_temp = mCitrine_stats(subset,:);
        
        replicate_list = unique(mCitrine_temp.replicate);
        duration_threshold_vs_WT = zeros(numel(replicate_list),1);
        
        % Calculate amplitude threshold vs WT
        y_diff_duration = zeros(numel(replicate_list),1);
        for kk = 1:numel(replicate_list)
            replicate = replicate_list(kk);
            
            mCitrine_temp_replicate = mCitrine_temp(mCitrine_temp.replicate==replicate,:);
            x = linspace(min(mCitrine_temp_replicate.mScarlet_AUC),max(mCitrine_temp_replicate.mScarlet_AUC),1000);
            y = interp1(mCitrine_temp_replicate.mScarlet_AUC,mCitrine_temp_replicate.mCitrine_max,x);
            
            [y_diff,x_50] = min(abs(y - duration_threshold_y_pMM0847));
            y_diff_duration(kk) = y_diff;
            duration_threshold_vs_WT(kk) = x(x_50);
            
        end
        duration_threshold_vs_WT_std = std(duration_threshold_vs_WT);
        duration_threshold_vs_WT = mean(duration_threshold_vs_WT);
        
        % Add value to table if distance to WT half max below threshold
        if sum(y_diff_amplitude<1)>1 && sum(y_diff_duration<1)>1
            promoter_thresholds_measured.amplitude_threshold_vs_WT(promoter_thresholds_measured.strain==strain & ...
                promoter_thresholds_measured.plasmid==plasmid) = amplitude_threshold_vs_WT;
            promoter_thresholds_measured.amplitude_threshold_vs_WT_std(promoter_thresholds_measured.strain==strain & ...
                promoter_thresholds_measured.plasmid==plasmid) = amplitude_threshold_vs_WT_std;
            
            promoter_thresholds_measured.duration_threshold_vs_WT(promoter_thresholds_measured.strain==strain & ...
                promoter_thresholds_measured.plasmid==plasmid) = duration_threshold_vs_WT;
            promoter_thresholds_measured.duration_threshold_vs_WT_std(promoter_thresholds_measured.strain==strain & ...
                promoter_thresholds_measured.plasmid==plasmid) = duration_threshold_vs_WT_std;
        end
        
    end
end

% Reset WT values
promoter_thresholds_measured.amplitude_threshold_vs_WT(promoter_thresholds_measured.plasmid=='pMM0847') = ...
    promoter_thresholds_measured.amplitude_threshold(promoter_thresholds_measured.plasmid=='pMM0847');
promoter_thresholds_measured.amplitude_threshold_vs_WT_std(promoter_thresholds_measured.plasmid=='pMM0847') = ...
    promoter_thresholds_measured.amplitude_threshold_std(promoter_thresholds_measured.plasmid=='pMM0847');
promoter_thresholds_measured.duration_threshold_vs_WT(promoter_thresholds_measured.plasmid=='pMM0847') = ...
    promoter_thresholds_measured.duration_threshold(promoter_thresholds_measured.plasmid=='pMM0847');
promoter_thresholds_measured.duration_threshold_vs_WT_std(promoter_thresholds_measured.plasmid=='pMM0847') = ...
    promoter_thresholds_measured.duration_threshold_std(promoter_thresholds_measured.plasmid=='pMM0847');

% Calculate mean of promoter thresholds across Msn2 mutants (excluding nonresponders)
% grp_vars = {'strain','reporter'};
% data_vars = {'amplitude_threshold','duration_threshold',...
%     'amplitude_threshold_vs_WT','duration_threshold_vs_WT',...
%     'slope_ratio','interpulse_decay'};
% nonresponders = (promoter_thresholds_measured.Msn2=='Msn2(WT|4E|T)' & ...
%     ismember(promoter_thresholds_measured.reporter,{'CTT1','SIP18','TKL2','ALD3','RTN2'}));
% promoter_thresholds_measured_mean = grpstats(promoter_thresholds_measured(~nonresponders,:),grp_vars,'nanmean','DataVars',data_vars);
% promoter_thresholds_measured_mean = clean_grpstats(promoter_thresholds_measured_mean);
% promoter_thresholds_measured_mean.Properties.VariableNames = regexprep(promoter_thresholds_measured_mean.Properties.VariableNames,'nanmean_','');

%% Calculate promoter thresholds from ideal pulses of nuclear Msn2

% Initialize variables
variables_list = {'strain','reporter','plasmid','Msn2','Msn2_tex','guess','k1','d1','k2','K','n','d2','k3','K_scale','fraction_active'};
subset = promoter_response_ideal.pulse_idx==1 & promoter_response_ideal.K_scale==1;
promoter_thresholds_ideal = unique(promoter_response_ideal(subset,variables_list));

promoter_thresholds_ideal.amp_thresh(:,1) = nan;
promoter_thresholds_ideal.amp_thresh_y(:,1) = nan;
promoter_thresholds_ideal.dur_thresh(:,1)= nan;
promoter_thresholds_ideal.dur_thresh_y(:,1)= nan;

%%%%%%%%%%%%%%%%%% Calculate thresholds per Msn2 mutant %%%%%%%%%%%%%%%%%%%
idx = 1;
for reporter_idx = 1:numel(reporter_list)
    reporter = reporter_list(reporter_idx);
    reporter = categorical(reporter);
    strain = unique(data_stats.strain(data_stats.reporter==reporter));
    disp(strain)
    
    % Plot fits over measurements
    for plasmid_idx = 1:numel(plasmid_list)
        
        plasmid = plasmid_list(plasmid_idx);
        Msn2 = unique(data_stats.Msn2(data_stats.plasmid==plasmid));
        Msn2_tex = unique(data_stats.Msn2_tex(data_stats.plasmid==plasmid));
        
        for guess = 1:n_guesses_plot
            
            % Calculate amplitude threshold
            promoter_params_temp = promoter_response_ideal(promoter_response_ideal.K_scale==1 & ...
                promoter_response_ideal.strain==strain & promoter_response_ideal.plasmid==plasmid & ...
                promoter_response_ideal.pulse_label=='amplitude' & promoter_response_ideal.A<=1 & ...
                promoter_response_ideal.guess==guess,:);
            
            x = linspace(min(promoter_params_temp.A),max(promoter_params_temp.A),1000);
            y = interp1(promoter_params_temp.A,promoter_params_temp.P_active_max,x);
            
            [~,amp_thresh] = min(abs(y - max(y)/2));
            amp_thresh_y = y(amp_thresh);
            amp_thresh = x(amp_thresh);
            amp_thresh = 100*amp_thresh;
            
            % Calculate duration threshold
            promoter_params_temp = promoter_response_ideal(promoter_response_ideal.K_scale==1 & ...
                promoter_response_ideal.strain==strain & promoter_response_ideal.plasmid==plasmid & ...
                promoter_response_ideal.pulse_label=='duration' & promoter_response_ideal.t1<=50 & ...
                promoter_response_ideal.guess==guess,:);
            
            x = linspace(min(promoter_params_temp.t1),max(promoter_params_temp.t1),1000);
            y = interp1(promoter_params_temp.t1,promoter_params_temp.P_active_max,x);
            
            [~,dur_thresh] = min(abs(y - max(y)/2));
            dur_thresh_y = y(dur_thresh);
            dur_thresh = x(dur_thresh);
            
            % Save thresholds
            promoter_thresholds_ideal.amp_thresh(promoter_thresholds_ideal.strain==strain & ...
                promoter_thresholds_ideal.plasmid==plasmid & promoter_thresholds_ideal.guess==guess) = amp_thresh;
            promoter_thresholds_ideal.amp_thresh_y(promoter_thresholds_ideal.strain==strain & ...
                promoter_thresholds_ideal.plasmid==plasmid & promoter_thresholds_ideal.guess==guess) = amp_thresh_y;
            promoter_thresholds_ideal.dur_thresh(promoter_thresholds_ideal.strain==strain & ...
                promoter_thresholds_ideal.plasmid==plasmid & promoter_thresholds_ideal.guess==guess) = dur_thresh;
            promoter_thresholds_ideal.dur_thresh_y(promoter_thresholds_ideal.strain==strain & ...
                promoter_thresholds_ideal.plasmid==plasmid & promoter_thresholds_ideal.guess==guess) = dur_thresh_y;
            idx = idx + 1;
        end
        
    end
    
end

%%%%%%%%%%%% Calculate amplitude/duration thresholds vs Msn2* %%%%%%%%%%%%%
promoter_thresholds_ideal.amp_thresh_vs_WT(:,1) = nan;
promoter_thresholds_ideal.amp_thresh_vs_WT_y(:,1) = nan;
promoter_thresholds_ideal.amp_thresh_vs_WT_rel(:,1) = nan;
promoter_thresholds_ideal.dur_thresh_vs_WT(:,1) = nan;
promoter_thresholds_ideal.dur_thresh_vs_WT_y(:,1) = nan;
promoter_thresholds_ideal.dur_thresh_vs_WT_rel(:,1) = nan;

idx = 1;
for reporter_idx = 1:numel(reporter_list)
    reporter = reporter_list(reporter_idx);
    reporter = categorical(reporter);
    strain = unique(data_stats.strain(data_stats.reporter==reporter));
    disp(strain)
    
    for plasmid_idx = 1:numel(plasmid_list)
        
        plasmid = plasmid_list(plasmid_idx);
        Msn2 = unique(data_stats.Msn2(data_stats.plasmid==plasmid));
        Msn2_tex = unique(data_stats.Msn2_tex(data_stats.plasmid==plasmid));
        
        % Calculate y value for half max P_active_k3 for amplitude conditions (mean across all parameter sets)
        promoter_params_temp = promoter_response_ideal(promoter_response_ideal.K_scale==1 & ...
            promoter_response_ideal.strain==strain & promoter_response_ideal.plasmid=='pMM0847' & ...
            promoter_response_ideal.pulse_label=='amplitude' & promoter_response_ideal.A<=1,:);
        promoter_params_temp = grpstats(promoter_params_temp,{'A'},'mean','DataVars','P_active_max_k3');
        promoter_params_temp = clean_grpstats(promoter_params_temp);
        promoter_params_temp.Properties.VariableNames(end) = {'P_active_max_k3'};
        
        x = linspace(min(promoter_params_temp.A),max(promoter_params_temp.A),1000);
        y = interp1(promoter_params_temp.A,promoter_params_temp.P_active_max_k3,x);
        
        [~,amp_thresh_WT] = min(abs(y - max(y)/2));
        amp_thresh_WT_y = y(amp_thresh_WT);
        amp_thresh_WT = x(amp_thresh_WT);
        amp_thresh_WT = 100*amp_thresh_WT;
        
        % Calculate y value for half max P_active_max_k3 for duration conditions (mean across all parameter sets)
        promoter_params_temp = promoter_response_ideal(promoter_response_ideal.K_scale==1 & ...
            promoter_response_ideal.strain==strain & promoter_response_ideal.plasmid=='pMM0847' & ...
            promoter_response_ideal.pulse_label=='duration' & promoter_response_ideal.t1<=50,:);
        promoter_params_temp = grpstats(promoter_params_temp,{'t1'},'mean','DataVars','P_active_max_k3');
        promoter_params_temp = clean_grpstats(promoter_params_temp);
        promoter_params_temp.Properties.VariableNames(end) = {'P_active_max_k3'};
             
        x = linspace(min(promoter_params_temp.t1),max(promoter_params_temp.t1),1000);        
        y = interp1(promoter_params_temp.t1,promoter_params_temp.P_active_max_k3,x);
               
        [~,dur_thresh_WT] = min(abs(y - max(y)/2));
        dur_thresh_WT_y = y(dur_thresh_WT);
        dur_thresh_WT = x(dur_thresh_WT);
        
        % Loop through all parameter sets and calculate amplitude threshold
        for guess = 1:n_guesses_plot
            
            % Calculate amplitude threshold vs WT
            promoter_params_temp = promoter_response_ideal(promoter_response_ideal.K_scale==1 & ...
                promoter_response_ideal.strain==strain & promoter_response_ideal.plasmid==plasmid & ...
                promoter_response_ideal.pulse_label=='amplitude' & promoter_response_ideal.A<=1 & ...
                promoter_response_ideal.guess==guess,:);
            
            x = linspace(min(promoter_params_temp.A),max(promoter_params_temp.A),1000);
            y = interp1(promoter_params_temp.A,promoter_params_temp.P_active_max_k3,x);
            
            [amp_thresh_vs_WT_diff,amp_thresh_vs_WT] = min(abs(y - amp_thresh_WT_y));
            amp_thresh_vs_WT_y = y(amp_thresh_vs_WT);
            amp_thresh_vs_WT = x(amp_thresh_vs_WT);
            amp_thresh_vs_WT = 100*amp_thresh_vs_WT;
            
            % Discard threshold if distance to WT half-max past threshold
            if amp_thresh_vs_WT_diff>1
                amp_thresh_vs_WT = nan;
            end
            
            amp_thresh_vs_WT_rel = amp_thresh_vs_WT/amp_thresh_WT;
            
            % Calculate duration threshold vs WT
            promoter_params_temp = promoter_response_ideal(promoter_response_ideal.K_scale==1 & ...
                promoter_response_ideal.strain==strain & promoter_response_ideal.plasmid==plasmid & ...
                promoter_response_ideal.pulse_label=='duration' & promoter_response_ideal.t1<=50 & ...
                promoter_response_ideal.guess==guess,:);
            
            % Predicted RTN2/TKL2 clearly overshoots measurements for 10
            % min pulses of Msn2(WT|4E|A) so forcing pulses <=10 min to
            % zero to avoid extreme underestimate of duration threshold            
            if ismember(reporter,{'RTN2','TKL2'})
                promoter_params_temp.P_active_max_k3(promoter_params_temp.t1<=10) = 0;
            end
            
            x = linspace(min(promoter_params_temp.t1),max(promoter_params_temp.t1),1000);
            y = interp1(promoter_params_temp.t1,promoter_params_temp.P_active_max_k3,x);
%             plot(x,y)
            
            [dur_thresh_vs_WT_diff,dur_thresh_vs_WT] = min(abs(y - dur_thresh_WT_y));
            dur_thresh_vs_WT_y = y(dur_thresh_vs_WT);
            dur_thresh_vs_WT = x(dur_thresh_vs_WT);
                        
            % Discard threshold if distance to WT half-max past threshold
            if dur_thresh_vs_WT_diff>1
                dur_thresh_vs_WT = nan;
            end
            
            dur_thresh_vs_WT_rel = dur_thresh_vs_WT/dur_thresh_WT;
            
            % Save amplitude and duration thresholds to table
            promoter_thresholds_ideal.amp_thresh_vs_WT(promoter_thresholds_ideal.strain==strain & ...
                promoter_thresholds_ideal.plasmid==plasmid & promoter_thresholds_ideal.guess==guess) = amp_thresh_vs_WT;
            promoter_thresholds_ideal.amp_thresh_vs_WT_y(promoter_thresholds_ideal.strain==strain & ...
                promoter_thresholds_ideal.plasmid==plasmid & promoter_thresholds_ideal.guess==guess) = amp_thresh_vs_WT_y;
            promoter_thresholds_ideal.amp_thresh_vs_WT_rel(promoter_thresholds_ideal.strain==strain & ...
                promoter_thresholds_ideal.plasmid==plasmid & promoter_thresholds_ideal.guess==guess) = amp_thresh_vs_WT_rel;
            promoter_thresholds_ideal.dur_thresh_vs_WT(promoter_thresholds_ideal.strain==strain & ...
                promoter_thresholds_ideal.plasmid==plasmid & promoter_thresholds_ideal.guess==guess) = dur_thresh_vs_WT;
            promoter_thresholds_ideal.dur_thresh_vs_WT_y(promoter_thresholds_ideal.strain==strain & ...
                promoter_thresholds_ideal.plasmid==plasmid & promoter_thresholds_ideal.guess==guess) = dur_thresh_vs_WT_y;
            promoter_thresholds_ideal.dur_thresh_vs_WT_rel(promoter_thresholds_ideal.strain==strain & ...
                promoter_thresholds_ideal.plasmid==plasmid & promoter_thresholds_ideal.guess==guess) = dur_thresh_vs_WT_rel;
            idx = idx + 1;
            
            
        end
    end
end

promoter_thresholds_ideal.amp_thresh_vs_WT(promoter_thresholds_ideal.plasmid=='pMM0847') = promoter_thresholds_ideal.amp_thresh(promoter_thresholds_ideal.plasmid=='pMM0847');
promoter_thresholds_ideal.dur_thresh_vs_WT(promoter_thresholds_ideal.plasmid=='pMM0847') = promoter_thresholds_ideal.dur_thresh(promoter_thresholds_ideal.plasmid=='pMM0847');

%%%%%%%%%%%%%%% Calculate slope ratio %%%%%%%%%%%%%%%%
% variables_list = {'strain','reporter','plasmid','Msn2','Msn2_tex','guess','params','K_scale','fraction_active'};

idx = 1;
for reporter_idx = 1:numel(reporter_list)
    reporter = reporter_list(reporter_idx);
    reporter = categorical(reporter);
    strain = unique(data_stats.strain(data_stats.reporter==reporter));
    disp(strain)
    
    % Plot fits over measurements
    for plasmid_idx = 1:numel(plasmid_list)
        
        plasmid = plasmid_list(plasmid_idx);
        Msn2 = unique(data_stats.Msn2(data_stats.plasmid==plasmid));
        Msn2_tex = unique(data_stats.Msn2_tex(data_stats.plasmid==plasmid));
        
        %%%%%%%%%%%%%%%%%%%%%% Calculate slope ratio %%%%%%%%%%%%%%%%%%%%%%

        % Continuous slope model
        promoter_params_temp = promoter_response_ideal(promoter_response_ideal.K_scale==1 & ...
            promoter_response_ideal.strain==strain & promoter_response_ideal.plasmid==plasmid & ...
            ismember(promoter_response_ideal.pulse_idx,[5,9,13]),:);
        
        mdl_continuous = fitglm(promoter_params_temp.mScarlet_AUC,promoter_params_temp.mCitrine_max);
        m_continuous = mdl_continuous.Coefficients.Estimate(2);
        m_continuous_se = mdl_continuous.Coefficients.SE(2);
        m_continuous_p = mdl_continuous.Coefficients.pValue(2);
        b_continuous = mdl_continuous.Coefficients.Estimate(1);
        
        x_mdl_continuous = (linspace(min(promoter_params_temp.mScarlet_AUC),max(promoter_params_temp.mScarlet_AUC),1000))';
        y_mdl_continuous = m_continuous.*x_mdl_continuous + b_continuous;
        
        % Pulsed slope model
        promoter_params_temp = promoter_response_ideal(promoter_response_ideal.K_scale==1 & ...
            promoter_response_ideal.strain==strain & promoter_response_ideal.plasmid==plasmid & ...
            ismember(promoter_response_ideal.pulse_idx,[35,36,37]),:);
        
        mdl_pulsed = fitglm(promoter_params_temp.mScarlet_AUC,promoter_params_temp.mCitrine_max);
        m_pulsed = mdl_pulsed.Coefficients.Estimate(2);
        m_pulsed_se = mdl_pulsed.Coefficients.SE(2);
        m_pulsed_p = mdl_pulsed.Coefficients.pValue(2);
        b_pulsed = mdl_pulsed.Coefficients.Estimate(1);
        
        x_mdl_pulsed = (linspace(min(promoter_params_temp.mScarlet_AUC),max(promoter_params_temp.mScarlet_AUC),1000))';
        y_mdl_pulsed = m_pulsed.*x_mdl_pulsed + b_pulsed;
        
        slope_ratio = m_pulsed./m_continuous;
        slope_ratio_se = (m_pulsed/m_continuous)* sqrt((m_pulsed_se/m_pulsed).^2 + (m_continuous_se/m_continuous).^2);
        
        % Interpulse model
        promoter_params_temp = promoter_response_ideal(promoter_response_ideal.K_scale==1 & ...
            promoter_response_ideal.strain==strain & promoter_response_ideal.plasmid==plasmid & ...
            ismember(promoter_response_ideal.pulse_idx,[9,36,38,39]),:);
        
        interpulse_function = @(a,x) a(1)*exp(a(2)*x); % a*exp(b*x)
        mdl_interpulse = fitnlm(promoter_params_temp.mScarlet_AUC,promoter_params_temp.mCitrine_max,interpulse_function,[1,0.1]);
        
        a_interpulse = mdl_interpulse.Coefficients.Estimate(1);
        b_interpulse = mdl_interpulse.Coefficients.Estimate(2);
        b_interpulse_se = mdl_interpulse.Coefficients.SE(2);
        b_interpulse_p = mdl_interpulse.Coefficients.pValue(2);
        
        x_mdl_interpulse = (linspace(min(promoter_params_temp.mScarlet_AUC),max(promoter_params_temp.mScarlet_AUC),1000))';
        y_mdl_interpulse = a_interpulse*exp(b_interpulse.*x_mdl_interpulse);
        
        % Save thresholds
        promoter_thresholds_ideal.slope_continuous(promoter_thresholds_ideal.strain==strain & ...
            promoter_thresholds_ideal.plasmid==plasmid) = m_continuous;
        promoter_thresholds_ideal.slope_continuous_se(promoter_thresholds_ideal.strain==strain & ...
            promoter_thresholds_ideal.plasmid==plasmid) = m_continuous_se;
        promoter_thresholds_ideal.slope_continuous_p(promoter_thresholds_ideal.strain==strain & ...
            promoter_thresholds_ideal.plasmid==plasmid) = m_continuous_p;
        
        promoter_thresholds_ideal.slope_pulsed(promoter_thresholds_ideal.strain==strain & ...
            promoter_thresholds_ideal.plasmid==plasmid) = m_pulsed;
        promoter_thresholds_ideal.slope_pulsed_se(promoter_thresholds_ideal.strain==strain & ...
            promoter_thresholds_ideal.plasmid==plasmid) = m_pulsed_se;
        promoter_thresholds_ideal.slope_pulsed_p(promoter_thresholds_ideal.strain==strain & ...
            promoter_thresholds_ideal.plasmid==plasmid) = m_pulsed_p;
        
        promoter_thresholds_ideal.slope_ratio(promoter_thresholds_ideal.strain==strain & ...
            promoter_thresholds_ideal.plasmid==plasmid) = slope_ratio;
        promoter_thresholds_ideal.slope_ratio_se(promoter_thresholds_ideal.strain==strain & ...
            promoter_thresholds_ideal.plasmid==plasmid) = slope_ratio_se;
        
        promoter_thresholds_ideal.interpulse_decay(promoter_thresholds_ideal.strain==strain & ...
            promoter_thresholds_ideal.plasmid==plasmid) = b_interpulse;
        promoter_thresholds_ideal.interpulse_decay_se(promoter_thresholds_ideal.strain==strain & ...
            promoter_thresholds_ideal.plasmid==plasmid) = b_interpulse_se;
        promoter_thresholds_ideal.interpulse_decay_p(promoter_thresholds_ideal.strain==strain & ...
            promoter_thresholds_ideal.plasmid==plasmid) = b_interpulse_p;        
    end    
end

%% Calculate AUC of ideal Msn2 trace (and ideal - measured)
t_Msn2_ideal = linspace(min(data_Msn2.time),max(data_Msn2.time),1000)';
t_measured = unique(data_Msn2.time);

idx = 1;
Msn2_ideal = cell(14,7);
for condition = 1:14
    Msn2_params = Msn2_params_list(Msn2_params_list.condition==condition,:);
    Msn2_localization_ideal = Msn2_CT(t_Msn2_ideal,Msn2_params);
    Msn2_AUC_ideal = trapz(t_Msn2_ideal,Msn2_localization_ideal);
    condition_display = unique(data_Msn2.condition_display(data_Msn2.condition==condition));    
    
    y_measured = data_Msn2.mScarlet_localization(data_Msn2.condition==condition);
    y_ideal = interp1(t_Msn2_ideal,Msn2_localization_ideal,t_measured);
    y_diff = y_ideal - y_measured;
    
    % Manually exlcude peaks where predicted and ideal offset by 1 frame
    exclude = y_diff<-0.025 | (t_measured<14 & y_diff>0.15);
    y_diff_exclude = y_diff;
    y_diff_exclude(exclude) = nan;
    
    Msn2_ideal{idx,1} = condition;
    Msn2_ideal{idx,2} = condition_display;
    Msn2_ideal{idx,3} = t_measured;
    Msn2_ideal{idx,4} = y_ideal;
    Msn2_ideal{idx,5} = y_diff;
    Msn2_ideal{idx,6} = y_diff_exclude;
    Msn2_ideal{idx,7} = Msn2_AUC_ideal;

    idx = idx + 1;
        
end

Msn2_ideal = cell2table(Msn2_ideal,'VariableNames',{'condition','condition_display','time','Msn2_localization_ideal',...
    'Msn2_localization_diff','Msn2_localization_exclude','Msn2_AUC_ideal'});

% Join ideal traces back into measurements table
data_Msn2_AUC = outerjoin(data_Msn2_AUC,Msn2_ideal(:,{'condition','Msn2_AUC_ideal'}), 'Type', 'Left', 'MergeKeys', true);

%% Cluster promoters by amp/dur thresholds and Kd

% Calc mean(log10(Kd)) for clustering
promoter_fits_mean = grpstats(promoter_fits(promoter_fits.rank<=100 & promoter_fits.Msn2=='Msn2(WT|4E|WT)',:),{'reporter'},@(x) mean(log10(x)),'DataVars','Kd');
promoter_fits_mean = clean_grpstats(promoter_fits_mean);
promoter_fits_mean.Properties.VariableNames(end) = {'log10_Kd'};

% Calculate mean thresholds per reporter and Msn2 mutant
promoter_thresholds_ideal_mean = grpstats(promoter_thresholds_ideal,{'strain','reporter','plasmid','Msn2','Msn2_tex','K_scale','fraction_active'},...
    'nanmean','DataVars',{'amp_thresh','dur_thresh','amp_thresh_vs_WT','amp_thresh_vs_WT_y','dur_thresh_vs_WT','dur_thresh_vs_WT_y','slope_ratio','interpulse_decay'});
promoter_thresholds_ideal_mean = clean_grpstats(promoter_thresholds_ideal_mean);
promoter_thresholds_ideal_mean.Properties.VariableNames = regexprep(promoter_thresholds_ideal_mean.Properties.VariableNames,'nanmean_','');

% Get oscillatory thresholds (slope ratio, interpulse decay) based on promoter response to measured pulses  of Msn2*
oscillatory_thresholds_WT = promoter_thresholds_measured(promoter_thresholds_measured.Msn2=='Msn2(WT|4E|WT)',{'reporter','slope_ratio','interpulse_decay'});

% Get amp/dur thresholds (slope ratio, interpulse decay) based on promoter response to ideal pulses of Msn2*
single_pulse_thresholds_WT = promoter_thresholds_ideal_mean(promoter_thresholds_ideal_mean.Msn2=='Msn2(WT|4E|WT)',{'reporter','amp_thresh','dur_thresh'});

% Merge oscillatory and amp/dur thresholds for Msn2* and cluster by them
promoter_thresholds_WT = outerjoin(oscillatory_thresholds_WT,single_pulse_thresholds_WT,'Type', 'Left', 'MergeKeys', true);
promoter_thresholds_WT = outerjoin(promoter_thresholds_WT,promoter_fits_mean,'Type', 'Left', 'MergeKeys', true);
promoter_thresholds_WT(promoter_thresholds_WT.reporter=='glpT',:)= [];
promoter_thresholds_WT.group = kmeans(zscore(promoter_thresholds_WT{:,{'amp_thresh','dur_thresh','log10_Kd'}}),3,'MaxIter',10000);

% Merge groups from clustering back into bigger table
promoter_thresholds_ideal_mean = outerjoin(promoter_thresholds_ideal_mean,promoter_thresholds_WT(:,{'reporter','group'}),'Type', 'Left', 'MergeKeys', true);

%% %%%%%%%%%%%%%%%% SAVE VARIABLES TO AVOID RECALCULATING %%%%%%%%%%%%%%%%%
clearvars -except data data_stats data_Msn2_AUC data_light_dose mCitrine_stats Msn2_ideal reporters_STRE_count...
    plot_colors promoter_fits promoter_fits_mean promoter_response_measured promoter_response_ideal ...
    promoter_thresholds_measured promoter_thresholds_measured_mean promoter_thresholds_ideal promoter_thresholds_ideal_mean

save('plot_figures_1_data','-v7.3')

return
%% %%%%%%%%%%%%%%%%%%%%%%%% PLOT MAIN FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 1B: plot HSP12 induction for Msn2 ± CLASP
close all;
Msn2_to_plot = {'Msn2','Msn2(WT|4E|WT)'};
% reporter_to_plot = 'HXK1';
reporter_to_plot = 'HSP12';

clc; clear g; figure('units','normalized','outerposition',[0 0 0.125 0.5]);
g = gramm('x',data_stats.time,'y',100*data_stats.mScarlet_localization,...
    'linestyle',data_stats.Msn2,'color',data_stats.Msn2,...
    'subset',data_stats.reporter==reporter_to_plot & ismember(data_stats.Msn2,Msn2_to_plot) & ismember(data_stats.condition,[1,4]));
g.facet_grid(data_stats.condition==1,[],'row_labels',false);
g.stat_summary('type','std','setylim',true);
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.set_names('x','time (min)','y',['Nuclear/cytoplasmic Msn2 (%)' newline 'mean ± std']);
g.set_names('x','','y','');
g.set_order_options('linestyle',-1,'row',0,'color',0);
g.set_text_options('interpreter','tex','base_size',12);
g.axe_property('XLim', [-10 120]);
% g.geom_polygon('x',{[0 20]},'color',[30 170 220]/255,'alpha',0.1);
% g.geom_polygon('x',{[0 0 20 20]},'y',{[-0.15 -0.075 -0.075 -0.15]},'color',[0 180 240]/255,'alpha',1);
g.no_legend();
g.draw();
% export_fig(fullfile(pwd,'localization_Msn2_and_Msn2(WT4EWT)'),'-png','-m4');

clc; clear g; figure('units','normalized','outerposition',[0 0 0.125 0.5]);
g = gramm('x',data_stats.time,'y',data_stats.mCitrine_cell,...
    'linestyle',data_stats.Msn2,...
    'subset',data_stats.reporter==reporter_to_plot & ismember(data_stats.Msn2,Msn2_to_plot) & ismember(data_stats.condition,[1,4]));
g.facet_grid(data_stats.condition==1,[],'row_labels',false);
g.stat_summary('type','std','setylim',true);
% g.set_color_options('map',[55 55 55]/255);
g.set_color_options('map',[230 170 0]/255);
% g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
% g.set_names('x','time (min)','y','HXK1-mCitrine')
g.set_names('x','','y','');
g.set_order_options('linestyle',-1,'row',0);
g.set_text_options('interpreter','tex','base_size',12);
g.axe_property('XLim', [-10 120]);
g.no_legend();
g.draw();
% export_fig(fullfile(pwd,'HSP12_induction_mCitrine_Msn2_and_Msn2(WT4EWT)'),'-png','-m4');

%% Figure 1D left panel: plot absolute nuclear mScarlet vs light dose
close all
Msn2_to_plot = {'Msn2','Msn2(WT|4E|WT)'};    

clc; clear g; figure('position',[100 100 400 300]);
g = gramm('x',data_light_dose.amplitudes,'y',(data_light_dose.mScarlet_nuclear),...
    'color',data_light_dose.Msn2_tex,'marker',data_light_dose.CLASP,...
    'subset',ismember(data_light_dose.Msn2,Msn2_to_plot)  & ...
    ismember(data_light_dose.bin,5:6) & data_light_dose.condition<=6);
g.stat_summary('type','std','geom',{'point','errorbar'},'setylim',true,'dodge',0.5,'width',2);
g.set_names('x','','y','','column','','row','','linestyle','','color','','marker','');
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.set_text_options('interpreter','tex','base_size',12);
g.set_line_options('style',{'-',':'});
g.set_point_options('markers',{'o','^'},'base_size',7);
g.axe_property('YLim',[7.5 57.5]);
g.no_legend();  
g.draw();
g.redraw(0.025);

set(findobj('color','g'),'Color',[85 85 85]/255);
% export_fig(fullfile(pwd,'Msn2_nuc_vs_light'),'-png','-m4');

%% Figure 1D right panel: plot absolute nuclear mScarlet vs light dose (PWM)
close all
Msn2_to_plot = {'Msn2','Msn2(WT|4E|WT)'};

clc; clear g; figure('position',[100 100 400 400]);
g = gramm('x',data_light_dose.PWM,'y',(data_light_dose.mScarlet_nuclear),...
    'color',data_light_dose.Msn2_tex,'marker',data_light_dose.CLASP,...
    'subset',ismember(data_light_dose.Msn2,Msn2_to_plot) & ...
    ismember(data_light_dose.bin,5:6) & data_light_dose.amplitudes==128);
g.stat_summary('type','std','geom',{'point','errorbar'},'setylim',true,'dodge',0.5,'width',2);
g.set_names('x','','y','','column','','row','','linestyle','','color','','marker','');
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.set_order_options('x',0,'color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.set_text_options('interpreter','tex','base_size',12);
g.set_line_options('style',{'-',':'});
g.set_point_options('markers',{'o','^'},'base_size',6);
g.axe_property('XTickLabelRotation',45,'YLim',[7.5 57.5]);
g.no_legend();
g.draw();
% export_fig(fullfile(pwd,'Msn2_nuc_vs_PWM'),'-png','-m4');

%% Figure 2A: heatmap of Msn2* light sweep experiments 

% Note: Msn2(WT|4E|A) = Msn2*


close all;
reporters_to_plot = {'glpT','RTN2','TKL2','SIP18','ALD3','CTT1','SIP18 D6','DCS2','SIP18 A4','DDR2','HXK1','HSP12'};
Msn2_list = categorical("Msn2(WT|4E|WT)");
cmap_red = zeros(255,3);
cmap_red(:,1) = 1:255;
cmap_red = cmap_red/255;

for mm = 1:numel(Msn2_list)
    
    clc; figure('units','normalized','outerposition',[0 0 0.75 0.3]);
    Msn2 = Msn2_list(mm);
    
    for rr = 1:(numel(reporters_to_plot)+1)
        subplot(1,numel(reporters_to_plot)+1,rr);
        
        if rr==1
            subset = ismember(data_stats.condition_display,1:14) & ismember(data_stats.Msn2,Msn2_list);
            
            data_heatmap = grpstats(data_stats(subset,:),{'condition_display','time','Msn2'},'nanmean','DataVars',{'mScarlet_localization'});
            data_heatmap = clean_grpstats(data_heatmap);
            data_heatmap.Properties.VariableNames(end) = {'mScarlet_localization'};
            data_heatmap = sortrows(data_heatmap,'condition_display');
            %             data_heatmap.mScarlet_nuclear = data_heatmap.mScarlet_nuclear./max(data_heatmap.mScarlet_nuclear);
            
            h = heatmap(data_heatmap,'time','condition_display','ColorVariable','mScarlet_localization','GridVisible','off','Colormap',cmap_red);
            h.Title = 'Msn2 (%)';
            
        else
            reporter = reporters_to_plot(rr-1);
            subset = ismember(data_stats.Msn2,Msn2) & data_stats.reporter==reporter & ismember(data_stats.condition_display,1:14);
            data_heatmap = grpstats(data_stats(subset,:),{'condition_display','time','Msn2'},'nanmean','DataVars',{'mCitrine_cell'});
            data_heatmap = clean_grpstats(data_heatmap);
            data_heatmap.Properties.VariableNames(end) = {'mCitrine_cell'};
            
            % Previously made mCitrine(mCitrine<1) = 1, but this messes up heat maps.
            % Subtracting 1 here for display purposes only
            data_heatmap.mCitrine_cell = (data_heatmap.mCitrine_cell - 1)./max(data_heatmap.mCitrine_cell);
            
            h = heatmap(data_heatmap(data_heatmap.Msn2==Msn2,:),'time','condition_display','ColorVariable','mCitrine_cell','GridVisible','off','Colormap',viridis);
            h.YDisplayLabels = nan(size(h.YDisplayData));
            h.YLabel = '';
            
            % Relabel negative control for display
            if ismember(reporter,'glpT')
                reporter = {'control'};
            end
            h.Title = reporter;
            
        end
        
        h.ColorLimits = [0 1];
        
        % Erase titles if not subplot(1,:)
        if mm~=1
            h.Title = '';
        end
        
        h.XDisplayLabels(1) = {'-10'};
        h.XDisplayLabels(end) = {'150'};
        h.XDisplayLabels(2:end-1) = {''};
        s = struct(h);
        s.XAxis.TickLabelRotation=0;
        
        % Add colorbar to only subplot(:,end)
        if rr~=(numel(reporters_to_plot)+1)
            colorbar('off');
        end
        
    end
    set(gcf,'color','w');
    
%     export_fig(fullfile(pwd,strcat(string(mm),'_light_sweep_heatmap_WT')),'-png','-m4');
end

%% Figure 2C: amplitude vs duration thresholds for Msn2* (ideal pulses)
close all
Msn2_to_plot = {'Msn2(WT|4E|WT)'};
nonresponders = (promoter_thresholds_ideal_mean.Msn2=='Msn2(WT|4E|T)' & ...
    ismember(promoter_thresholds_ideal_mean.reporter,{'CTT1','SIP18','TKL2','ALD3','RTN2'}));

clc; clear g; figure('position',[100 100 700 700]);
g = gramm('x',promoter_thresholds_ideal_mean.dur_thresh,'y',promoter_thresholds_ideal_mean.amp_thresh,...
    'label',promoter_thresholds_ideal_mean.reporter,...
    'subset',promoter_thresholds_ideal_mean.K_scale==1 & promoter_thresholds_ideal_mean.reporter~='glpT' & ...
    ismember(promoter_thresholds_ideal_mean.Msn2,Msn2_to_plot) & ~nonresponders);
g.facet_wrap(promoter_thresholds_ideal_mean.Msn2_tex,'ncols',3);
g.geom_point();
g.set_names('x','','y','','column','','color','');
g.set_order_options('column',0);
g.set_color_options('map',[55 55 55]/255);
g.set_text_options('base_size',14,'interpreter','tex');
g.set_point_options('base_size',8);
g.no_legend();
g.draw();

% Add promoter labels
g.update('x',promoter_thresholds_ideal_mean.dur_thresh + 0.75,'label',promoter_thresholds_ideal_mean.reporter);
g.geom_label('FontSize',14,'FontWeight', 'Bold');
g.no_legend();
g.draw();
% export_fig(fullfile(pwd,'amp_vs_dur_thresh'),'-png','-m4');

%% Figure 2D top panel: Histogram of predicted Kd for RTN2 and HSP12
close all
clc; clear g; figure('position',[0 0 550 275]);
reporters_to_plot = {'RTN2','HSP12'};

g = gramm('x',log10(promoter_fits.Kd),...
    'color',promoter_fits.rank<=100,...
    'subset',promoter_fits.Msn2=='Msn2(WT|4E|WT)' & log10(promoter_fits.Kd)<10 & ...
    ismember(promoter_fits.reporter,reporters_to_plot));
g.facet_grid([],promoter_fits.reporter,'row_labels',false),
g.stat_bin('geom','overlaid_bar','normalization','pdf','nbins',10,'fill','transparent')
g.set_names('x','','y','','color','','column','');
g.set_order_options('column',reporters_to_plot);
g.set_color_options('map',[125 125 125; 0 0 0;]/255);
g.set_text_options('interpreter','tex','base_size',14);
g.no_legend();
g.axe_property('XTickLabelRotation',45);
g.draw()
g.redraw(0.025);
% export_fig(fullfile(pwd,'Kd_RTN2_HSP12'),'-png','-m4');

% Test if K^n of top 0.1% of RTN2 parameters significantly different than bottom 99.9%
RTN2_top = log10(promoter_fits.Kd(promoter_fits.Msn2=='Msn2(WT|4E|WT)' & promoter_fits.reporter=='RTN2' & promoter_fits.rank<=100));
RTN2_rest = log10(promoter_fits.Kd(promoter_fits.Msn2=='Msn2(WT|4E|WT)' & promoter_fits.reporter=='RTN2' & promoter_fits.rank>100));
[h p1] = kstest2(RTN2_top,RTN2_rest)

% Test if K^n of top 0.1% of HSP12 parameters significantly different than bottom 99.9%
HSP12_top = log10(promoter_fits.Kd(promoter_fits.Msn2=='Msn2(WT|4E|WT)' & promoter_fits.reporter=='HSP12' & promoter_fits.rank<=100));
HSP12_rest = log10(promoter_fits.Kd(promoter_fits.Msn2=='Msn2(WT|4E|WT)' & promoter_fits.reporter=='HSP12' & promoter_fits.rank>100));
[h p2] = kstest2(HSP12_top,HSP12_rest)

%% Figure 2D bottom panel: Predicted Kd for all reporters
close all
reporters_to_plot = {'RTN2','TKL2','SIP18','ALD3','CTT1','SIP18 D6','DCS2','SIP18 A4','DDR2','HXK1','HSP12'};
Msn2_to_plot = {'Msn2(WT|4E|WT)'};
n_guesses_plot = 100;
nonresponders = (promoter_fits.Msn2=='Msn2(WT|4E|T)' & ...
    ismember(promoter_fits.reporter,{'CTT1','SIP18','TKL2','ALD3','RTN2'}));

clc; clear g; figure('position',[100 100 600 350]);
g = gramm('x',promoter_fits.reporter,'y',log10(promoter_fits.K.^promoter_fits.n),...
    'color',promoter_fits.Msn2_tex,...
    'subset',ismember(promoter_fits.Msn2,Msn2_to_plot) & ...
    ismember(promoter_fits.reporter,reporters_to_plot) & ...
    promoter_fits.rank<=n_guesses_plot & ~nonresponders);
g.stat_summary('type','ci','geom',{'bar','black_errorbar'},'setylim',true,'width',0.7,'dodge',0.75);
g.set_names('x','','y','','color','','column','');
g.axe_property('XTickLabelRotation',45);
g.set_text_options('interpreter','tex','base_size',14,'big_title_scaling',1);
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.set_order_options('x',reporters_to_plot);
g.set_color_options('map',[125 125 125]/255);
g.axe_property('YLim',[0 4.2]);
g.no_legend();
g.draw();
% export_fig(fullfile(pwd,'Kd_all_reporters'),'-png','-m4');

%% Figure 3C: RTN2, DCS2, and HSP12 time courses w/ K scaling (ideal pulses, WT only)

%%%%%%%%%%% CHANGED TO PROMOTER_RESPONSE_MEASURED %%%%%%%%%%%%%%%

close all
reporters_to_plot = {'RTN2','DCS2','HSP12'};
Msn2_to_plot = {'Msn2(WT|4E|A)','Msn2(WT|4E|WT)','Msn2(WT|4E|T)'};

% clc; clear g; figure('position',[100 100 1600 500]);
clc; clear g; figure('position',[100 100 700 225]);
g = gramm('x',promoter_response_measured.time,'y',promoter_response_measured.mCitrine,...
    'color',promoter_response_measured.K_scale,...
    'subset',ismember(promoter_response_measured.reporter,reporters_to_plot) & ...
    ismember(promoter_response_measured.Msn2,'Msn2(WT|4E|WT)') & ...
    promoter_response_measured.condition==9);
g.facet_wrap(promoter_response_measured.reporter,'ncols',6,'scale','independent');
g.stat_summary('type','ci','setylim',true);
g.set_names('x','','y','','column','','color','');
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.set_order_options('column',reporters_to_plot,'color',[0.5 1 2]);
g.set_text_options('base_size',11,'interpreter','tex');
g.no_legend();
g.draw();
% g.facet_axes_handles(1).YLim = [-5 50]
% export_fig(fullfile(pwd,'mCitrine_K_scaling'),'-png','-m4');

% Calculate stats to quantify these differences in paper
% promoter_response_ideal_mean = promoter_response_measured;
% subset = ismember(promoter_response_ideal_mean.reporter,{'RTN2','DCS2','HSP12'}) & ...
%     promoter_response_ideal_mean.Msn2=='Msn2(WT|4E|WT)' & ...
%     promoter_response_ideal_mean.pulse_idx==29;
% promoter_response_ideal_mean = promoter_response_ideal_mean(subset,:);
% grp_vars = {'reporter','Msn2','K_scale','pulse_idx'};
% promoter_response_ideal_mean = grpstats(promoter_response_ideal_mean,grp_vars,'mean','DataVars','mCitrine_max');
% promoter_response_ideal_mean = clean_grpstats(promoter_response_ideal_mean);
% promoter_response_ideal_mean.Properties.VariableNames(end) = {'mCitrine_max'};
% 
% promoter_response_ideal_mean_K_scale_1 = promoter_response_ideal_mean(promoter_response_ideal_mean.K_scale==1,{'reporter','pulse_idx','mCitrine_max'});
% promoter_response_ideal_mean_K_scale_1.Properties.VariableNames(end) = {'mCitrine_max_K_scale_1'};
% promoter_response_ideal_mean = outerjoin(promoter_response_ideal_mean,promoter_response_ideal_mean_K_scale_1, 'Type', 'Left', 'MergeKeys', true);
% 
% promoter_response_ideal_mean.mCitrine_max_vs_K_scale_1 = promoter_response_ideal_mean.mCitrine_max./promoter_response_ideal_mean.mCitrine_max_K_scale_1;

%% Figure 3D: plot max predicted reporter expression for select localization patterns

promoter_response_mean = promoter_response_measured(promoter_response_measured.plasmid=='pMM0847',:);
promoter_response_mean = grpstats(promoter_response_mean,{'strain','reporter','plasmid','condition','K_scale'},'mean','DataVars','mCitrine_max');
promoter_response_mean = clean_grpstats(promoter_response_mean);
promoter_response_mean.Properties.VariableNames('mean_mCitrine_max') = {'mCitrine_max'};

promoter_response_mean_Kscale1 = promoter_response_mean(promoter_response_mean.K_scale==1,{'strain','reporter','condition','mCitrine_max'});
promoter_response_mean_Kscale1.Properties.VariableNames('mCitrine_max') = {'mCitrine_max_Kscale1'};
promoter_response_mean = outerjoin(promoter_response_mean,promoter_response_mean_Kscale1,'Type', 'Left', 'MergeKeys', true);
promoter_response_mean.condition_label(:,1) = categorical("NA");

promoter_response_mean.condition_label(promoter_response_mean.condition==9) = categorical("continuous");
promoter_response_mean.condition_label(promoter_response_mean.condition==10) = categorical("pulsed");
promoter_response_mean.condition_label(promoter_response_mean.condition==3) = categorical("low");
promoter_response_mean.condition_label(promoter_response_mean.condition==2) = categorical("short");


reporters_to_plot = {'RTN2','DCS2','HSP12'};

close all
clc; clear g; figure('position',[100 100 700 225]);
g = gramm('x',(promoter_response_mean.condition_label),'y',(promoter_response_mean.mCitrine_max),...
    'color',(promoter_response_mean.K_scale),'subset',ismember(promoter_response_mean.condition,[2,3,9,10]));
g.facet_grid([],promoter_response_mean.reporter,'scale','independent','row_labels',false);
g.set_names('x','','y','mCitrine','column','','color','');
g.geom_bar('width',0.5,'LineWidth',0.25)
g.set_color_options('map',[7 182	75; 0 160 240; 240 135 50]/255);
g.set_order_options('x',{'short','low','continuous','pulsed'},'column',reporters_to_plot,'color',0);
g.draw();
% export_fig(fullfile(pwd,'fig_2C_new_absolute'),'-png','-m4');

clc; clear g; figure('position',[100 100 700 225]);
g = gramm('x',(promoter_response_mean.condition_label),'y',(promoter_response_mean.mCitrine_max./promoter_response_mean.mCitrine_max_Kscale1),...
    'color',(promoter_response_mean.K_scale),'subset',ismember(promoter_response_mean.condition,[2,3,9,10]));
g.facet_grid([],promoter_response_mean.reporter,'scale','independent','row_labels',false);
g.set_names('x','','y','relative mCitrine','column','','color','');
g.geom_bar('width',0.5,'LineWidth',0.25);
g.set_color_options('map',[7 182	75; 0 160 240; 240 135 50]/255);
g.set_order_options('x',{'short','low','continuous','pulsed'},'column',reporters_to_plot,'color',0);
% g.axe_property('YLim',[0 10],'YTickLabel','')
g.axe_property('YLim',[0 2.5],'YTickLabel','')
g.draw();
% Note: add page break manually by compositing plots with different YLims
% g.facet_axes_handles(1).YLim = [0 10];
% g.facet_axes_handles(3).YLim = [0.7 1.25];
% export_fig(fullfile(pwd,'fig_2C_new_relative'),'-png','-m4');
% export_fig(fullfile(pwd,'fig_2C_new_relative_break'),'-png','-m4');

%% Figure 3E: expression vs time for 50 min 100% pulse Msn2(A)*, Msn2*, and Msn2(T)* (measurements)

close all
reporters_to_plot = {'RTN2','DCS2','HSP12'};
Msn2_to_plot = {'Msn2(WT|4E|A)','Msn2(WT|4E|WT)','Msn2(WT|4E|T)'};

clc; clear g; figure('position',[100 100 700 225]);
g = gramm('x',data_stats.time,'y',data_stats.mCitrine_cell,...
    'color',data_stats.Msn2_tex,...
    'subset',ismember(data_stats.reporter,reporters_to_plot) & ...
    ismember(data_stats.Msn2,Msn2_to_plot) & ...
    data_stats.condition==9);
g.facet_wrap(data_stats.reporter,'ncols',3,'scale','independent');
g.stat_summary('type','std','setylim',true);
g.set_names('x','','y','','column','','color','');
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.set_order_options('column',reporters_to_plot,'color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.set_text_options('base_size',11,'interpreter','tex');
g.no_legend();
g.draw();
% export_fig(fullfile(pwd,'mCitrine_mutants_measured'),'-png','-m4');

%% Figure 4: heatmaps for Msn2(A)*, Msn2*, and Msn2(T)* light sweep experiments

close all;
reporters_to_plot = {'glpT','RTN2','TKL2','SIP18','ALD3','CTT1','SIP18 D6','DCS2','SIP18 A4','DDR2','HXK1','HSP12'};
Msn2_list = categorical({'Msn2(WT|4E|A)','Msn2(WT|4E|WT)','Msn2(WT|4E|T)'});
% cmap_red = [(1:255)' zeros(255,2)]/255;
% cmap_red(cmap_red<10) = 0;
cmap_red = zeros(255,3);
cmap_red(:,1) = 1:255;
cmap_red = cmap_red/255;

for mm = 1:numel(Msn2_list)
    
    clc; figure('units','normalized','outerposition',[0 0 0.75 0.3]);
    Msn2 = Msn2_list(mm);
    
    for rr = 1:(numel(reporters_to_plot)+1)
        subplot(1,numel(reporters_to_plot)+1,rr);
        
        if rr==1
            subset = ismember(data_stats.condition,1:14) & ismember(data_stats.Msn2,Msn2_list);
            
            data_heatmap = grpstats(data_stats(subset,:),{'condition_display','time','Msn2'},'nanmean','DataVars',{'mScarlet_localization'});
            data_heatmap = clean_grpstats(data_heatmap);
            data_heatmap.Properties.VariableNames(end) = {'mScarlet_localization'};
            %             data_heatmap.mScarlet_nuclear = data_heatmap.mScarlet_nuclear./max(data_heatmap.mScarlet_nuclear);
            
            h = heatmap(data_heatmap,'time','condition_display','ColorVariable','mScarlet_localization','GridVisible','off','Colormap',cmap_red);
            h.Title = 'Msn2 (%)';
            
        else
            reporter = reporters_to_plot(rr-1);
            subset = data_stats.reporter==reporter & ismember(data_stats.condition,1:14);
            data_heatmap = grpstats(data_stats(subset,:),{'condition_display','time','Msn2'},'nanmean','DataVars',{'mCitrine_cell'});
            data_heatmap = clean_grpstats(data_heatmap);
            data_heatmap.Properties.VariableNames(end) = {'mCitrine_cell'};
            
            % Previously made mCitrine(mCitrine<1) = 1, but this messes up heat maps.
            % Subtracting 1 here for display purposes only
            data_heatmap.mCitrine_cell = (data_heatmap.mCitrine_cell - 1)./max(data_heatmap.mCitrine_cell);
            
            h = heatmap(data_heatmap(data_heatmap.Msn2==Msn2,:),'time','condition_display','ColorVariable','mCitrine_cell','GridVisible','off','Colormap',viridis);
            h.YDisplayLabels = nan(size(h.YDisplayData));
            h.YLabel = '';
            
            % Relabel negative control for display
            if ismember(reporter,'glpT')
                reporter = {'control'};
            end
            
            h.Title = reporter;
            
        end
        
        h.ColorLimits = [0 1];
        
        % Erase titles if not subplot(1,:)
        if mm~=1
            h.Title = '';
        end
        
        h.XDisplayLabels(1) = {'-10'};
        h.XDisplayLabels(end) = {'150'};
        h.XDisplayLabels(2:end-1) = {''};
        s = struct(h);
        s.XAxis.TickLabelRotation=0;
        
        % Add colorbar to only subplot(:,end)
        if rr~=(numel(reporters_to_plot)+1)
            colorbar('off');
        end
        
    end
    set(gcf,'color','w');
%     export_fig(fullfile(pwd,strcat(string(mm),'_light_sweep_heatmap')),'-png','-m4');
end

%% Figures 5A - 5B: plot mCitrine vs amplitude/duration/slope ratio conditions (measured, all mutants)
close all
Msn2_to_plot = {'Msn2(WT|4E|A)','Msn2(WT|4E|WT)','Msn2(WT|4E|T)'};
reporters_to_plot = {'RTN2','HSP12'};

% Plot mCitrine vs Msn2 AUC for amplitude conditions
amplitude_conditions = 1:2:9;
clc; clear g; figure('position',[100 100 250 400]);
g = gramm('x',(mCitrine_stats.mScarlet_AUC),'y',mCitrine_stats.mCitrine_max,...
    'color',(mCitrine_stats.Msn2_tex),...
    'subset',ismember(mCitrine_stats.condition,amplitude_conditions) ...
    & ismember(mCitrine_stats.Msn2,Msn2_to_plot) & ismember(mCitrine_stats.reporter,reporters_to_plot));
g.facet_grid(mCitrine_stats.reporter,[],'row_labels',false,'scale','independent');
g.stat_summary('type','std','geom',{'point','errorbar'},'setylim',true,'width',2,'dodge',0.75);
g.set_names('x','','y','','column','','row','','linestyle','','color','','marker','');
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.set_order_options('row',reporters_to_plot,'color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.set_text_options('base_size',11,'interpreter','tex','title_scaling',1);
g.set_point_options('base_size',6,'markers',{'o','^'});
g.axe_property('XLim',[-5 45]);
g.no_legend();
g.draw()
% export_fig(fullfile(pwd,'mCitrine_vs_amplitude'),'-png','-m4');

% Plot duration conditions
duration_conditions = [1,2:2:8,9];
clc; clear g; figure('position',[100 100 250 400]);
g = gramm('x',mCitrine_stats.mScarlet_AUC,'y',mCitrine_stats.mCitrine_max,...
    'color',mCitrine_stats.Msn2_tex,...
    'subset',ismember(mCitrine_stats.condition,duration_conditions) & ...
    ismember(mCitrine_stats.Msn2,Msn2_to_plot) & ismember(mCitrine_stats.reporter,reporters_to_plot));
g.facet_grid(mCitrine_stats.reporter,[],'row_labels',false,'scale','independent');
g.stat_summary('type','std','geom',{'point','errorbar'},'setylim',true,'width',4,'dodge',0.75);
g.set_names('x','Msn2 AUC','y','max mCitrine','column','','row','','linestyle','','color','','marker','increasing duration)');
g.set_names('x','','y','','column','','row','','linestyle','','color','','marker','');
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.set_order_options('row',reporters_to_plot,'color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.set_text_options('base_size',10,'interpreter','tex','title_scaling',1);
g.set_point_options('base_size',6,'markers',{'o','^'});
g.axe_property('XLim',[-5 45]);
g.no_legend();
g.draw()
% export_fig(fullfile(pwd,'mCitrine_vs_duration'),'-png','-m4');

% Plot slope ratio for continuous vs pulsed conditions
continuous_conditions = [1,2,6,9];
clc; clear g; figure('position',[0 0 600 425]);
g = gramm('x',mCitrine_stats.mScarlet_AUC,'y',mCitrine_stats.mCitrine_max,...
    'color',mCitrine_stats.Msn2_tex,...
    'subset',ismember(mCitrine_stats.condition,continuous_conditions) & ...
    ismember(mCitrine_stats.Msn2,Msn2_to_plot) & ismember(mCitrine_stats.reporter,reporters_to_plot));
g.facet_grid(mCitrine_stats.reporter,mCitrine_stats.Msn2,'scale','free_y','row_labels',false,'column_labels',false);
g.stat_summary('type','std','geom',{'point','errorbar'},'setylim',true,'width',0.75);
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.set_order_options('linestyle',0,'column',Msn2_to_plot,'row',reporters_to_plot,'color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.set_text_options('base_size',10,'interpreter','tex','title_scaling',1);
% g.set_names('x','Msn2 AUC','y','max mCitrine','color','','row','','column','','marker','continuous(1) vs pulsed(0)','linestyle','continuous(1) vs pulsed(0)');
g.set_names('x','','y','','color','','row','','column','','marker','continuous(1) vs pulsed(0)','linestyle','continuous(1) vs pulsed(0)');
g.set_point_options('markers',{'o'},'base_size',6);
g.set_line_options('styles',{'-'},'base_size',1.5);
g.stat_fit('fun',@(m,b,x) m*x + b,'intopt','functional');
g.no_legend();
g.draw();

for ii = 1:size(g.results.stat_fit,1)
    g.results.stat_fit(ii).area_handle.FaceAlpha = 0.06;
end

pulsed_conditions = [1,10:12];
g.update('subset',ismember(mCitrine_stats.condition,pulsed_conditions) & ...
    ismember(mCitrine_stats.Msn2,Msn2_to_plot) & ismember(mCitrine_stats.reporter,reporters_to_plot));
g.stat_summary('type','std','geom',{'point','errorbar'},'width',0.75);
g.set_color_options('map',0.85*plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.stat_fit('fun',@(m,b,x) m*x + b,'intopt','functional');
g.no_legend();
g.set_point_options('markers',{'^'},'base_size',5.5);
g.set_line_options('styles',{':'},'base_size',1.5);
g.draw();
for ii = 1:size(g.results.stat_fit,1)
    g.results.stat_fit(ii).area_handle.FaceAlpha = 0.06;
end
% export_fig(fullfile(pwd,'mCitrine_pulsed_vs_continuous'),'-png','-m4');

%% Figure 5C: plot slope ratio w/ variance
reporters_to_plot = {'glpT','RTN2','TKL2','SIP18','ALD3','CTT1','SIP18 D6','DCS2','SIP18 A4','DDR2','HXK1','HSP12'};
Msn2_to_plot = {'Msn2(WT|4E|A)','Msn2(WT|4E|WT)','Msn2(WT|4E|T)'};
nonresponders = (promoter_thresholds_measured.Msn2=='Msn2(WT|4E|T)' & ...
    ismember(promoter_thresholds_measured.reporter,{'CTT1','SIP18','TKL2','ALD3','RTN2'}));

close all
clc; clear g; figure('position',[100 100 500 250]);
g = gramm('x',promoter_thresholds_measured.reporter,'y',promoter_thresholds_measured.slope_ratio,...
    'ymin',promoter_thresholds_measured.slope_ratio - promoter_thresholds_measured.slope_ratio_se,...
    'ymax',promoter_thresholds_measured.slope_ratio + promoter_thresholds_measured.slope_ratio_se,...
    'color',promoter_thresholds_measured.Msn2_tex,...
    'subset',promoter_thresholds_measured.reporter~='glpT' & ismember(promoter_thresholds_measured.Msn2,Msn2_to_plot) & ~nonresponders);
g.geom_interval('geom',{'bar','black_errorbar'},'width',0.7,'dodge',0.75);
g.set_names('x','','y','','column','','color','');
g.set_order_options('x',reporters_to_plot,'color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.set_text_options('interpreter','tex','base_size',11);
g.no_legend();
g.axe_property('XTickLabelRotation',45);
g.draw();
g.redraw(0.05);
% export_fig(fullfile(pwd,'slope_ratios'),'-png','-m4');

%% Figure 5D: plot mCitrine histograms

exclude = data.reporter=='ALD3' & data.Msn2=='Msn2(WT|4E|A)' & data.replicate=='2';
Msn2_to_plot = {'Msn2(WT|4E|A)','Msn2(WT|4E|WT)','Msn2(WT|4E|T)'};
reporter_list = {'RNT2','HSP12'}

for r = 1:numel(reporter_list)
    reporter = char(reporter_list(r));
    experiments_to_plot = unique(data_stats.experiment(data_stats.condition==1 & data_stats.reporter==reporter));
    
    subset = ismember(data.experiment,experiments_to_plot) & ismember(data.Msn2,'Msn2(WT|4E|T)') & ...
        data.condition==1 & data.time>115 & data.time<125;
    threshold_dark = prctile(data.mCitrine_cell(subset,:),99);

    subset = ismember(data.experiment,experiments_to_plot) & ismember(data.Msn2,'Msn2(WT|4E|T)') & ...
        data.condition==6 & data.time>115 & data.time<125;
    threshold_continuous = sum(data.mCitrine_cell(subset,:)>threshold_dark)./size(data(subset,:),1);

    subset = ismember(data.experiment,experiments_to_plot) & ismember(data.Msn2,'Msn2(WT|4E|T)') & ...
        data.condition==11 & data.time>115 & data.time<125;
    threshold_pulsed = sum(data.mCitrine_cell(subset,:)>threshold_dark)./size(data(subset,:),1);

    close all;
    clc; clear g; figure('position',[100 100 550 150]);
    g = gramm('x',log10(data.mCitrine_cell),...
        'color',cellstr(data.Msn2_tex),...
        'subset',ismember(data.experiment,experiments_to_plot) & ismember(data.Msn2,Msn2_to_plot) & ...
        ismember(data.condition,[1,6,11]) & (data.time>115 & data.time<125));
    g.facet_wrap((data.condition),'ncols',3,'column_labels',0);
    g.stat_bin('nbins',12,'geom','stairs','normalization','pdf');
    g.set_names('x','','y','','column','','row','','linestyle','','color','');    
    g.set_order_options('color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
    g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
    g.set_text_options('base_size',10,'interpreter','tex','title_scaling',1);
    g.no_legend();
    g.geom_vline('xintercept',log10(threshold_dark))
    g.draw();
    disp(strcat("T_continous: ", num2str(100*threshold_continuous),"% > threshold"))
    disp(strcat("T_pulsed: ", num2str(100*threshold_pulsed),"% > threshold"))
%     export_fig(fullfile(pwd,[reporter '_pulsed_vs_continuous_hist_stairs']),'-png','-m4');

end
%% %%%%%%%%%%%%%%%%%%%% PLOT SUPPLEMENTARY FIGURES %%%%%%%%%%%%%%%%%%%%%%%%
%% Figure S1E: plot basal mCitrine histograms
close all
Msn2_to_plot = {'Msn2','Msn2(WT|4E|WT)'};

data_basal = data(ismember(data.condition,1:15) & ...
    data.time<0 & data.mScarlet_nuclear>0 & data.mScarlet_nuclear<=65,:);

clc; clear g; figure('position',[100 100 500 300]);
g = gramm('x',data_basal.mScarlet_nuclear,'color',data_basal.Msn2_tex,'subset',data_basal.condition<=14);
g.stat_bin('geom','line','normalization','pdf','nbins',16);
g.set_names('x','Nuclear Msn2 in dark','y','Density','column','','row','','linestyle','','color','','marker','');
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.set_order_options('color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.set_text_options('interpreter','tex','base_size',11);
g.no_legend();
g.draw();

clc; clear g; figure('position',[100 100 500 300]);
g = gramm('x',data_basal.mScarlet_nuclear,'color',data_basal.Msn2_tex,'subset',data_basal.condition==15);
g.stat_bin('geom','line','normalization','pdf','nbins',16);
g.set_names('x','Nuclear Msn2 in dark','y','Density','column','','row','','linestyle','','color','','marker','');
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.set_order_options('color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.set_text_options('interpreter','tex','base_size',11);
g.no_legend();
g.draw();

numel(unique(data_basal.experiment(data_basal.plasmid=='pMM0845'))) % # experiments Msn2-CLASP
numel(unique(data_basal.experiment(data_basal.plasmid=='pMM0847'))) % # experiments Msn2*-CLASP

%% Figure S1F: plot endpoint mCitrine measurements for Msn2-CLASP and Msn2*-CLASP
close all

% Add simplified labels
data_stats.Msn2_label_simple(:,1) = "NA";
data_stats.Msn2_label_simple(data_stats.plasmid=='pMM0845') = "Msn2-CLASP";
data_stats.Msn2_label_simple(data_stats.plasmid=='pMM0847') = "Msn2*-CLASP";
data_stats.Msn2_label_simple = categorical(data_stats.Msn2_label_simple);

reporters_to_plot = {'glpT','TKL2','ALD3','SIP18 D6','DCS2','HXK1','RTN2','SIP18','CTT1','SIP18 A4','DDR2','HSP12'};
Msn2_to_plot = {'Msn2','Msn2(WT|4E|WT)'};

clc; clear g; figure('position',[100 100 1000 450]);
g = gramm('x',data_stats.Msn2_label_simple,'y',data_stats.mCitrine_cell,...
    'color',data_stats.Msn2_tex,...
    'subset',ismember(data_stats.Msn2,Msn2_to_plot) & data_stats.time>120 & ...
    data_stats.condition==1);
g.facet_wrap(data_stats.reporter,'ncols',6,'scale','independent');
g.stat_summary('type','std','geom',{'bar','black_errorbar'},'width',1.5);
g.set_color_options('map',0.85*plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.set_order_options('x',-1,'column',reporters_to_plot,'color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.set_text_options('interpreter','tex');
g.set_names('x','','y','mCitrine (AU)','column','','color','');
g.no_legend();
g.draw();
g.redraw(0.06);

% Adjust YLim and export
for ii = 1:12
    y_lim = g.facet_axes_handles(ii).YLim;
    if y_lim(2)<10
        y_lim(2) = 10;
    end
    g.facet_axes_handles(ii).YLim = [-0.25 y_lim(2)];
end
% export_fig(fullfile(pwd,'endpoint_mCitrine_condition1_Msn2_vs_Msn2_WT4EWT'),'-png','-m4');


%% Figure S2A: plot light sweep timecourses
close all;
reporters_to_plot = {'glpT','TKL2','ALD3','SIP18 D6','DCS2','HXK1','RTN2','SIP18','CTT1','SIP18 A4','DDR2','HSP12'};

%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot light programs %%%%%%%%%%%%%%%%%%%%%%%%%%%
%Import example light program
led_program = load(fullfile(pwd,'example_led_program.mat'));
led_program = led_program.experiment;
led_program = get_led_patterns(led_program,160*60,true);
led_program.well = categorical(led_program.well);
led_program = led_program(ismember(led_program.column,{'01','02'}),:);
led_program.time = cellfun(@(x) (x/60) - 10,led_program.time,'UniformOutput',false);
well_conditions = unique(data_stats(:,{'well','condition_display'}),'rows');
led_program = outerjoin(led_program,well_conditions, 'Type', 'Left', 'MergeKeys', true);

% Plot light program
clc; clear g; figure('position',[100 100 2400 125]);
g = gramm('x',led_program.time(:),'y',led_program.intensity(:),'group',led_program.condition_display);
g.facet_wrap(led_program.condition_display,'ncols',14,'column_labels',false,'scale','fixed');
g.geom_line();
g.set_names('x','time (min.)','y','light (AU)','column','','row','','linestyle','','color','');
g.set_color_options('map',[86 128 202]/255);
% g.set_title(strcat('LED program'));
g.no_legend();
g.draw()
g.redraw(0.025);
% export_fig(fullfile(pwd,'light_programs'),'-png','-m4');

%%%%%%%%%%%%%% Plot composite of Msn2 localization per mutant %%%%%%%%%%%%%
Msn2_to_plot = {'Msn2(WT|4E|A)','Msn2(WT|4E|WT)','Msn2(WT|4E|T)'};
clc; clear g; figure('position',[100 100 2400 125]);
g = gramm('x',data_stats.time,'y',100*data_stats.mScarlet_localization,...
    'color',data_stats.Msn2_tex,...
    'subset',data_stats.condition_display<=14 & ismember(data_stats.Msn2,Msn2_to_plot) & ...
    ~ismember(data_stats.reporter,{'ALD3','CTT1'}));
g.facet_wrap(data_stats.condition_display,'ncols',14,'column_labels',false);
g.stat_summary('type','std','setylim',true);
g.set_names('x','time (min)','y','Msn2 amplitude (%)','color','','column','');
g.set_text_options('base_size',10,'interpreter','tex','title_scaling',1);
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.set_order_options('row',reporters_to_plot,'color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.no_legend();
g.draw();
g.redraw(0.025);
% export_fig(fullfile(pwd,'Msn2_loc_all_mutants'),'-png','-m4');

%%%%%%%%%%%%%%%%%%%% Plot fits over mCitrine timecourses %%%%%%%%%%%%%%%%%%
close all;
conditions_to_plot = 1:14;
Msn2_to_plot ={'Msn2(WT|4E|A)','Msn2(WT|4E|WT)','Msn2(WT|4E|T)'};

reporter_list = unique(data_stats.reporter);
for ii = 1:numel(reporter_list)
    
    close all;
    reporter = reporter_list(ii);
    
    % Plot measured mCitrine
    clc; clear g; figure('position',[100 100 2400 125]);
    g = gramm('x',data_stats.time,'y',data_stats.mCitrine_cell,...
        'color',data_stats.Msn2_tex,...
        'subset',ismember(data_stats.reporter,reporter) & ismember(data_stats.Msn2,Msn2_to_plot) &...
        ismember(data_stats.condition,conditions_to_plot));
    g.facet_wrap(data_stats.condition_display,'ncols',14,'column_labels',false);
    g.stat_summary('type','std');
    g.set_point_options('base_size',2);
    g.set_names('x','time (min.)','y',strcat('p',string(reporter),'-mCitrine'),'column','','row','','color','');
    g.set_text_options('interpreter','tex');
    g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
    g.set_order_options('color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
%     g.set_title(string(reporter));
    g.set_line_options('base_size',1);
    g.no_legend();
    g.draw();
    
    % Plot simulated mCitrine
    g.update('x',promoter_response_measured.time,'y',promoter_response_measured.mCitrine,...
        'color',promoter_response_measured.Msn2_tex,'subset',ismember(promoter_response_measured.reporter,reporter) & ...
        ismember(promoter_response_measured.Msn2,Msn2_to_plot) & ismember(promoter_response_measured.condition,conditions_to_plot) & ...
        promoter_response_measured.K_scale==1);
    g.facet_wrap(promoter_response_measured.condition_display,'ncols',14);
    g.stat_summary('type','std');
    g.set_color_options('map',0.5*plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
    g.set_order_options('color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
    g.set_line_options('styles',{':'});
    g.no_legend();
    g.draw();
    g.redraw(0.025);
%     export_fig(fullfile(pwd,strcat(string(reporter),'_mCitrine_vs_time_with_fits')),'-png','-m4');
end
%% Figure S2C: plot reporter induction for full light CLASP vs dCLASP

close all;
Msn2_to_plot = {'Msn2(WT|4E|A)','Msn2(WT|4E|WT)','Msn2(WT|4E|T)'};
reporters_to_plot = {'glpT','RTN2','TKL2','SIP18','ALD3','CTT1','SIP18 D6','DCS2','SIP18 A4','DDR2','HXK1','HSP12'};

% Exclude outlier replicate
exclude = (mCitrine_stats.reporter=='ALD3' & mCitrine_stats.Msn2=='Msn2(WT|4E|A)' & mCitrine_stats.replicate=='2');

clc; clear g; figure('units','normalized','outerposition',[0.1 0.1 0.45 0.5]);
g = gramm('x',categorical(mCitrine_stats.CLASP),'y',mCitrine_stats.mCitrine_max,...
    'color',mCitrine_stats.Msn2_tex,...
    'marker',mCitrine_stats.CLASP,...
    'subset',ismember(mCitrine_stats.Msn2,Msn2_to_plot) & ismember(mCitrine_stats.condition,[9,15]) & ~exclude);
g.facet_wrap(mCitrine_stats.reporter,'ncols',4,'scale','independent');
g.stat_summary('type','std','geom',{'bar','black_errorbar'},'setylim',true);
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.set_order_options('column',reporters_to_plot,'x',0,'color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.set_names('x','','y','max mCitrine (AU)','column','')
g.set_text_options('interpreter','tex','base_size',8);
g.axe_property('XLabelScale',0.1);
g.no_legend();
g.draw();
g.redraw(0.05);
g.facet_axes_handles(1).YLim = [-5 50];

%% Figure S3A: Plot top 100 (0.1%) of each parameter for all promoters
close all

reporters_to_plot = {'RTN2','TKL2','SIP18','ALD3','CTT1','SIP18 D6','DCS2','SIP18 A4','DDR2','HXK1','HSP12'};
promoter_fits.Kd = promoter_fits.K.^promoter_fits.n;
params_to_plot = {'k1','d1','k2','K','n','Kd','d2','k3'};
Msn2_to_plot = {'Msn2(WT|4E|A)','Msn2(WT|4E|WT)','Msn2(WT|4E|T)'};
Msn2_to_plot = {'Msn2(WT|4E|WT)'}

nonresponders = (promoter_fits.Msn2=='Msn2(WT|4E|T)' & ...
    ismember(promoter_fits.reporter,{'CTT1','SIP18','TKL2','ALD3','RTN2'}));

clc; clear g; figure('position',[0 0 1800 250]);
for ii = 1:numel(params_to_plot)
    parameter_to_plot = string(params_to_plot(ii));
    
    if strcmp(parameter_to_plot,'n')
        y = promoter_fits.(parameter_to_plot);
        y_label = parameter_to_plot;
    elseif strcmp(parameter_to_plot,'Kd')
        y = log10(promoter_fits.(parameter_to_plot));
        y_label = strcat("log_{10}(K^{n})")
    else
        y = log10(promoter_fits.(parameter_to_plot));
        y_label = strcat("log_{10}(",parameter_to_plot,")")
    end
    
    g(1,ii) = gramm('x',promoter_fits.reporter,'y',y,...
        'color',promoter_fits.Msn2_tex,...
        'subset',ismember(promoter_fits.Msn2,Msn2_to_plot) & promoter_fits.rank<=100 & ...
        ismember(promoter_fits.reporter,reporters_to_plot) & ~nonresponders);
    
    g(1,ii).stat_summary('type','ci','geom',{'point','errorbar'},'setylim',true,'width',1.5);
    g(1,ii).set_names('x','','y',y_label,'color','','column','');
    g(1,ii).set_order_options('x',reporters_to_plot);
    g(1,ii).set_color_options('map',[55 55 55]/255);
    g(1,ii).set_text_options('interpreter','tex','base_size',10);
    g(1,ii).axe_property('XTickLabelRotation',60);
    g(1,ii).no_legend();
end
g.draw()
% export_fig(fullfile(pwd,'top_100_params'),'-png','-m4');
%% Figure S3B: measured vs ideal Msn2 vs time (composite)

%%%%%%%%%%%% Plot measured vs ideal localization time courses %%%%%%%%%%%%%
close all
% Msn2_ideal = outerjoin(Msn2_ideal,condition_reorder, 'Type', 'Left', 'MergeKeys', true);

% Plot ideal Msn2 localization
clc; clear g; figure('units','normalized','outerposition',[0 0 0.7 0.5]);
g = gramm('x',Msn2_ideal.time,'y',cellfun(@(x) 100.*x, Msn2_ideal.Msn2_localization_ideal,'Un',0));
g.facet_wrap(Msn2_ideal.condition_display,'ncols',7);
g.geom_line();
g.set_names('x','time (min.)','y','Rel. nuclear Msn2 (%)','column','');
g.set_color_options('map',[125 125 125]/255);
g.set_text_options('Interpreter','tex');
g.set_line_options('styles',{':'},'base_size',1)
g.set_title('Measured vs ideal Msn2 localization');

% Plot measured localization
g.update('x',data_stats.time,'y',100*data_stats.mScarlet_localization,...
    'subset',data_stats.condition<=14);
g.facet_wrap(data_stats.condition_display,'ncols',7);
g.stat_summary('type','std','setylim',true);
g.set_color_options('map',[55 55 55]/255);
g.set_line_options('styles',{'-'},'base_size',1.5)
g.draw();
% export_fig(fullfile(pwd,'Msn2_localization_measured_vs_ideal'),'-png','-m4');

numel(unique(data_stats.experiment(data_stats.plasmid=='pMM1079'))) % # experiments Msn2(T)*-CLASP

%%%%%%%%%%%%%%% Plot AUC for measured vs ideal localization %%%%%%%%%%%%%%%

% Plot measured vs ideal Msn2 AUC
clc; clear g; figure('position',[100 100 300 300]);
g = gramm('x',data_Msn2_AUC.Msn2_AUC_ideal,'y',data_Msn2_AUC.mScarlet_AUC,...
    'subset',data_Msn2_AUC.condition_display<=14);
g.geom_point();
g.geom_abline();
% g.set_names('x','Msn2 AUC (ideal)','y','Msn2 AUC (measured)','column','');
g.set_names('x','','y','','column','');
g.set_color_options('map',[55 55 55]/255);
g.set_point_options('markers',{'o','^'});
g.axe_property('XLim',[0 50],'YLim',[0 50]);
g.set_text_options('Interpreter','tex','base_size',12);
g.draw();

g.update('x',data_Msn2_AUC.Msn2_AUC_ideal + 1,'label',data_Msn2_AUC.condition_display);
g.geom_label();
g.draw();
% export_fig(fullfile(pwd,'Msn2_AUC_measured_vs_ideal'),'-png','-m4');

%% Figure S3E: plot amp_thresh_vs_WT and dur_thresh_vs_WT calculation example

close all

Msn2_to_plot = {'Msn2(WT|4E|A)','Msn2(WT|4E|WT)','Msn2(WT|4E|T)'};
reporter = 'SIP18 A4';

% Plot amplitude threshold
clc; clear g; figure('position',[100 100 400 300]);
g = gramm('x',100*promoter_response_ideal.A,'y',promoter_response_ideal.P_active_max_k3,...
    'color',promoter_response_ideal.Msn2_tex,'group',promoter_response_ideal.Msn2,...
    'subset',ismember(promoter_response_ideal.pulse_idx,19:29) & promoter_response_ideal.reporter==reporter & promoter_response_ideal.K_scale==1);
g.facet_wrap(promoter_response_ideal.reporter,'ncols',4,'scale','independent');
g.set_names('x','amplitude (%)','y','k_{3}P_{on}','column','','color','');
g.set_order_options('color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.stat_summary('type','std');
g.set_point_options('base_size',2);
g.set_text_options('interpreter','tex','base_size',11);
g.no_legend();
g.draw();
g.redraw(0.05);

% Plot activation timescale
clc; clear g; figure('position',[100 100 400 300]);
g = gramm('x',promoter_response_ideal.t1,'y',promoter_response_ideal.P_active_max_k3,...
    'color',promoter_response_ideal.Msn2_tex,...
    'subset',ismember(promoter_response_ideal.pulse_idx,1:13) & promoter_response_ideal.reporter==reporter & promoter_response_ideal.K_scale==1);
g.facet_wrap(promoter_response_ideal.reporter,'ncols',4,'scale','independent');
g.set_names('x','duration (min)','y','k_{3}P_{on}','column','','color','');
g.set_order_options('column',1:11,'color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.stat_summary('type','std');
g.set_point_options('base_size',2);
g.set_text_options('interpreter','tex','base_size',11);
g.no_legend();
g.draw();
g.redraw(0.05);

%% Figure S3F: generic plot of transition rate as a function of Msn2 with varying K and n

close all; clc

% Define Hill function used to model promoter transition
Hill = @(TF,k1,n,K) (k1.*TF.^n)./(K.^n + TF.^n);
x = 0:0.01:10;

% Plot transition rate vs K
figure; hold on
plot(x,Hill(x,1,1,2),'m')
plot(x,Hill(x,1,1,1),'k')
plot(x,Hill(x,1,1,0.5),'r')
plot(x,Hill(x,1,1,0.25),'g')
plot(x,0.5*ones(length(x)),'k--') % hline at half-max
legend('K = 2','K = 1','K = 0.5','K = 0.25','Location','northwest')
title('k1 = 1, n = 1')
ylim([0 1.1])
xlim([0 1.1])
xlabel('[Msn2]')
ylabel('transition rate')

% Plot transition rate vs K
figure; hold on
plot(x,Hill(x,1,0.5,0.5),'m')
plot(x,Hill(x,1,1,0.5),'k')
plot(x,Hill(x,1,2,0.5),'r')
plot(x,Hill(x,1,4,0.5),'g')
plot(x,0.5*ones(length(x)),'k--') % hline at half-max
legend('n = 0.5','n = 1','n = 2','n = 4', 'Location','northwest')
title('k1 = 1, K = 0.5')
ylim([0 1.1])
xlim([0 1.1])
xlabel('[Msn2]')
ylabel('transition rate')

%% Figure S3G: plot predicted Kd vs # STREs

% load('reporters_STRE_count')
% Add STRE counts to promoter fits tables
promoter_fits = outerjoin(promoter_fits,reporters_STRE_count, 'Type', 'Left', 'MergeKeys', true);
promoter_fits_mean = outerjoin(promoter_fits_mean,reporters_STRE_count, 'Type', 'Left', 'MergeKeys', true);

close all
reporters_to_plot = {'RTN2','TKL2','SIP18','ALD3','CTT1','SIP18 D6','DCS2','SIP18 A4','DDR2','HXK1','HSP12'};
Msn2_to_plot = {'Msn2(WT|4E|WT)'};
n_guesses_plot = 100;
nonresponders = (promoter_fits.Msn2=='Msn2(WT|4E|T)' & ...
    ismember(promoter_fits.reporter,{'CTT1','SIP18','TKL2','ALD3','RTN2'}));

clc; clear g; figure('position',[100 100 700 400]);
g = gramm('x',promoter_fits.STRE_count_1000,'y',log10(promoter_fits.Kd),...
    'group',promoter_fits.reporter,...
    'subset',ismember(promoter_fits.Msn2,Msn2_to_plot) & ...
    ismember(promoter_fits.reporter,reporters_to_plot) & ...
    promoter_fits.rank<=n_guesses_plot & ~nonresponders);
g.stat_summary('type','ci','geom',{'point','errorbar'});
g.set_names('x','# STREs','y',['log_{10}(K^{n})' newline '(mean ± ci)'],'color','','column','');
g.set_names('x','','y','','color','','column','');
g.set_text_options('interpreter','tex','base_size',11,'big_title_scaling',1 );
g.set_color_options('map',[55 55 55]/255);
g.axe_property('YLim',[0 4.2]);
g.no_legend();
g.draw();

% Add reporter labels
g.update('x',promoter_fits_mean.STRE_count_1000 + 0.075,'y',promoter_fits_mean.log10_Kd,...
    'label',promoter_fits_mean.reporter,'group',promoter_fits_mean.reporter,...
    'subset',ismember(promoter_fits_mean.reporter,reporters_to_plot));
g.geom_label();
g.set_names('x','# STREs','y','mean ± 95% CI of log_{10}(K^{n})','color','','column','');
g.draw();
% export_fig(fullfile(pwd,'Kd_vs_STREs_1000bp_upstream_all_reporters'),'-png','-m4');


%% Figure S4G: plot amplitude and duration thresholds vs WT w/ variance (ideal pulses) ****
reporters_to_plot = {'glpT','RTN2','TKL2','SIP18','ALD3','CTT1','SIP18 D6','DCS2','SIP18 A4','DDR2','HXK1','HSP12'};
Msn2_to_plot = {'Msn2(WT|4E|A)','Msn2(WT|4E|WT)','Msn2(WT|4E|T)'};
nonresponders = (promoter_thresholds_ideal.Msn2=='Msn2(WT|4E|T)' & ...
    ismember(promoter_thresholds_ideal.reporter,{'CTT1','SIP18','TKL2','ALD3','RTN2'}));

close all
clc; clear g; figure('position',[100 100 800 350]);
g = gramm('x',promoter_thresholds_ideal.reporter,'y',promoter_thresholds_ideal.amp_thresh_vs_WT,...
    'color',promoter_thresholds_ideal.Msn2_tex,...
    'subset',promoter_thresholds_ideal.reporter~='glpT' & ~nonresponders);
g.stat_summary('type','ci','geom',{'bar','black_errorbar'},'width',0.7,'dodge',0.75);
g.set_names('x','','y','amplitude threshold (%)','column','','color','');
g.set_order_options('x',reporters_to_plot,'color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.set_text_options('interpreter','tex','base_size',11,'big_title_scaling',1);
g.no_legend();
g.axe_property('XTickLabelRotation',45,'YLim',[0 105]);
g.draw();
g.redraw(0.05);

clc; clear g; figure('position',[100 100 800 350]);
g = gramm('x',promoter_thresholds_ideal.reporter,'y',(promoter_thresholds_ideal.dur_thresh_vs_WT),...
    'color',promoter_thresholds_ideal.Msn2_tex,...
    'subset',promoter_thresholds_ideal.reporter~='glpT' & ~nonresponders);
g.stat_summary('type','ci','geom',{'bar','black_errorbar'},'width',0.7,'dodge',0.75);
g.set_names('x','','y','activation timescale (min)','column','','color','');
g.set_order_options('x',reporters_to_plot,'color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.set_text_options('interpreter','tex','base_size',11,'big_title_scaling',1);
g.no_legend();
g.axe_property('XTickLabelRotation',45);
g.draw();
g.redraw(0.05);

%% Figure S5A: compare the SIPs
close all
reporters_to_plot = {'SIP18','SIP18 D6','SIP18 A4'};
Msn2_to_plot = {'Msn2(WT|4E|A)','Msn2(WT|4E|WT)','Msn2(WT|4E|T)'};

clc; clear g; figure('position',[100 100 800 400]);
g = gramm('x',data_stats.time,'y',data_stats.mCitrine_cell,...
    'color',data_stats.Msn2_tex,...
    'subset',ismember(data_stats.reporter,reporters_to_plot) & ...
    ismember(data_stats.Msn2,Msn2_to_plot) & ...
    ismember(data_stats.condition,[6,11]));
g.facet_grid(data_stats.condition,data_stats.reporter,'scale','independent','row_labels',false);
g.stat_summary('type','std','setylim',true);
g.set_names('x','','y','' ,'column','','color','');
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.set_order_options('column',reporters_to_plot,'color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.set_text_options('base_size',11,'interpreter','tex');
g.no_legend();
g.draw();

% Update YLims and export
g.facet_axes_handles(1,1).YLim = [-5 100];
g.facet_axes_handles(2,1).YLim = [-5 100];
g.facet_axes_handles(1,2).YLim = [-15 250];
g.facet_axes_handles(2,2).YLim = [-15 250];
g.facet_axes_handles(1,3).YLim = [-25 625];
g.facet_axes_handles(2,3).YLim = [-25 625];
% export_fig(fullfile(pwd,'compare_SIPs'),'-png','-m4');

%% Figure S5B: plot mCitrine noise 

close all

subset = data_stats.time>120;
grp_vars = {'strain','plasmid','condition','replicate'};
mCitrine_noise = grpstats(data_stats(subset,:),grp_vars,'nanmean','DataVars','mCitrine_noise');
mCitrine_noise = clean_grpstats(mCitrine_noise);
mCitrine_noise.Properties.VariableNames(end) = {'mCitrine_noise'};
mCitrine_stats = outerjoin(mCitrine_stats,mCitrine_noise, 'Type', 'Left', 'MergeKeys', true);

Msn2_to_plot = {'Msn2(WT|4E|A)','Msn2(WT|4E|WT)','Msn2(WT|4E|T)'};
nonresponders = (mCitrine_stats.Msn2=='Msn2(WT|4E|T)' & ...
    ismember(mCitrine_stats.reporter,{'CTT1','SIP18','TKL2','ALD3','RTN2'}));
exclude = mCitrine_stats.reporter=='HXK1' & mCitrine_stats.Msn2=='Msn2(WT|4E|T)' & ...
    mCitrine_stats.condition==11 & mCitrine_stats.replicate=='1';
reporters_to_plot = {'RTN2','TKL2','SIP18','ALD3','CTT1','SIP18 D6','DCS2','SIP18 A4','DDR2','HXK1','HSP12'};

% Plot noise for copntinuous conditions
clc; clear g; figure('position',[100 100 500 250]);
g = gramm('x',(mCitrine_stats.reporter),'y',(mCitrine_stats.mCitrine_noise),...
    'color',mCitrine_stats.Msn2_tex,'group',mCitrine_stats.replicate,...
    'subset',ismember(mCitrine_stats.reporter,reporters_to_plot) & ismember(mCitrine_stats.Msn2,Msn2_to_plot) & ...
    ismember(mCitrine_stats.condition,6) & ~exclude);
g.stat_summary('type','std','geom',{'point'},'width',0.7,'dodge',0.75,'setylim',true);
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.set_order_options('x',reporters_to_plot,'color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.set_names('x','','y','','color','','column','','row','');
g.set_point_options('markers',{'o'},'base_size',6);
g.axe_property('XTickLabelRotation',45,'YLim',[-0.5 10]);
g.set_text_options('interpreter','tex','base_size',11);
g.no_legend();
g.draw()

% Plot noise for pulsed conditions
g.update('x',(mCitrine_stats.reporter),'y',(mCitrine_stats.mCitrine_noise),...
    'color',mCitrine_stats.Msn2_tex,'group',mCitrine_stats.replicate,...
    'subset',ismember(mCitrine_stats.reporter,reporters_to_plot) & ismember(mCitrine_stats.Msn2,Msn2_to_plot) & ...
    ismember(mCitrine_stats.condition,11) & ~exclude);
g.stat_summary('type','std','geom',{'point'},'width',0.7,'dodge',0.75,'setylim',true);
g.set_color_options('map',0.85*plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.set_order_options('x',reporters_to_plot,'color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.set_point_options('markers',{'^'},'base_size',6);
g.no_legend();
g.draw()
% export_fig(fullfile(pwd,'mCitrine_noise'),'-png','-m4');

%% Figure S5C: max mCitrine vs noise
grp_vars = {'strain','reporter','plasmid','Msn2','Msn2_tex','condition'};
mCitrine_max_mean = grpstats(mCitrine_stats,grp_vars,'mean','DataVars',{'mCitrine_max','mCitrine_noise'});
mCitrine_max_mean = clean_grpstats(mCitrine_max_mean);
mCitrine_max_mean.Properties.VariableNames(end-1:end) = {'mCitrine_max','mCitrine_noise'};
nonresponders = (mCitrine_max_mean.Msn2=='Msn2(WT|4E|T)' & ...
    ismember(mCitrine_max_mean.reporter,{'CTT1','SIP18','TKL2','ALD3','RTN2'}));

close all
clc; clear g; figure('position',[100 100 475 450]);
g = gramm('x',log10(mCitrine_max_mean.mCitrine_max),'y',log10(mCitrine_max_mean.mCitrine_noise),...
    'color',mCitrine_max_mean.Msn2_tex,'label',mCitrine_max_mean.reporter,...
    'subset',ismember(mCitrine_max_mean.reporter,reporters_to_plot) & ismember(mCitrine_max_mean.Msn2,Msn2_to_plot) & ...
    mCitrine_max_mean.condition==9 & ~nonresponders);
g.geom_point();
g.geom_label();
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.set_order_options('x',reporters_to_plot,'color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.set_names('x','log_{10}(max mCitrine)','y','log_{10}(noise)','color','','column','','row','');
g.set_text_options('interpreter','tex','base_size',11);
g.no_legend();
g.draw();

%% %%%%%%%%%%%%%%%%%%%%%%%% MISCELLANEOUS PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Simulate mRNA expression to determine sampling time for FISH experiments

close all
reporters_to_plot = {'RTN2','HSP12'};
Msn2_to_plot = {'Msn2(WT|4E|A)','Msn2(WT|4E|WT)','Msn2(WT|4E|T)'};

clc; clear g; figure('position',[100 100 1200 400]);
g = gramm('x',promoter_response_measured.time,'y',promoter_response_measured.mRNA,...
    'color',promoter_response_measured.Msn2_tex,...
    'subset',ismember(promoter_response_measured.reporter,reporters_to_plot) & ...
    ismember(promoter_response_measured.Msn2,Msn2_to_plot) & ...
    promoter_response_measured.K_scale==1 & promoter_response_measured.fraction_active==1 & ...
    promoter_response_measured.guess<=10 & promoter_response_measured.condition==9);
g.facet_grid([],promoter_response_measured.reporter,'scale','independent')
g.stat_summary();
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.set_order_options('x',reporters_to_plot,'color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.set_text_options('Interpreter','tex');
g.set_names('x','time (min)','y','predicted mRNA','color','','column','');
g.draw();

%% Plot mCitrine levels normalized to Msn2*
close all

mCitrine_stats_mean = grpstats(mCitrine_stats,{'strain','plasmid','reporter','Msn2_tex','condition'},'mean','DataVars','mCitrine_max');
mCitrine_stats_mean.Properties.VariableNames('mean_mCitrine_max') = {'mCitrine_max'};
% 
mCitrine_stats_mean_WT = mCitrine_stats_mean(mCitrine_stats_mean.plasmid=='pMM0847',:);
mCitrine_stats_mean_WT.Properties.VariableNames('mCitrine_max') = {'mCitrine_max_WT'};
mCitrine_stats_mean_WT.plasmid = [];
mCitrine_stats_mean_WT.Msn2_tex = [];

mCitrine_stats_mean = outerjoin(mCitrine_stats_mean,mCitrine_stats_mean_WT, 'Type', 'Left', 'MergeKeys', true);
mCitrine_stats = outerjoin(mCitrine_stats,mCitrine_stats_mean_WT, 'Type', 'Left', 'MergeKeys', true);

Msn2_to_plot = {'Msn2(WT|4E|A)','Msn2(WT|4E|WT)','Msn2(WT|4E|T)'};


close all
clc; clear g; figure('position',[100 100 900 300]);
g = gramm('x',categorical(mCitrine_stats.condition),'y',(mCitrine_stats.mCitrine_max./mCitrine_stats.mCitrine_max_WT),...
    'color',mCitrine_stats.Msn2_tex,'subset',ismember(mCitrine_stats.condition,[3,9]) & ...
    ismember(mCitrine_stats.plasmid,{'pMM0846','pMM0847','pMM1079'}) & ismember(mCitrine_stats.reporter,{'HSP12'}));
g.facet_wrap(mCitrine_stats.reporter,'ncols',6)
g.stat_summary('type','std','geom',{'bar','black_errorbar'});
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.set_order_options('color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.set_text_options('Interpreter','tex');
g.set_names('x','','y','max mCitrine vs Msn2*','color','')
% g.axe_property('XLim',[1 11]);
g.draw();

g.update();
g.geom_point('dodge',0.6)  
g.set_color_options('map',[155 155 155; 155 155 155; 155 155 155]/255)
g.draw()
