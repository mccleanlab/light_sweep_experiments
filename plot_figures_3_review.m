close all; clearvars; clc;

% This script simulates the expression of hypothetical promoters while
% scaling the parameters of the gene expression model one-by-one. 
%% %%%%%%%%%%%%%%%%%%%% IMPORT DATA AND RUN SIMULATIONS %%%%%%%%%%%%%%%%%%%
%% Select input/output folders

parent_folder  = pwd;
if contains(parent_folder,'\')
    parent_folder_idx   = strfind(parent_folder,'\');
else
    parent_folder_idx   = strfind(parent_folder,'/');
end
parent_folder = parent_folder(1:parent_folder_idx(end)-1);

%% Import data and parameters

% Import LHS promoter parameters
model_folder = fullfile(parent_folder,'promoter_model');
model_solutions_folder = fullfile(model_folder,'output_1_Hill_1_Msn2_100k_d2');
load(fullfile(model_solutions_folder,'promoter_params_LHS.mat'));

% Load measurements
load(fullfile(parent_folder,'light_sweep_experiments','mCitrine_stats.mat'))
load(fullfile(parent_folder,'light_sweep_experiments','data_stats.mat'))
load(fullfile(parent_folder,'light_sweep_experiments','data_Msn2.mat'))
load(fullfile(parent_folder,'light_sweep_experiments','data_Msn2_AUC.mat'))
load(fullfile(parent_folder,'promoter_analysis','reporter_promoters.mat'))

% Process
data_stats = data_stats(data_stats.condition<=14,:);
data_stats = data_stats(ismember(data_stats.plasmid,{'pMM0846','pMM0847','pMM1079'}),:);
data_stats.mCitrine_cell = data_stats.mCitrine_cell - data_stats.mCitrine_cell_basal;

% Get plot colors and Msn2_tex labels
opts = detectImportOptions(fullfile(parent_folder,'plot_settings.xlsx'),'Sheet','Msn2_tex');
opts.VariableTypes = {'categorical','categorical','double','double','double'};
plot_colors = readtable(fullfile(parent_folder,'plot_settings.xlsx'),opts);
plot_colors.RGB = [plot_colors.R, plot_colors.G, plot_colors.B];
plot_colors = plot_colors(:,{'Msn2','Msn2_tex','RGB'});
plot_colors.RGB(plot_colors.Msn2=='Msn2',:) = [70, 115, 190];
Msn2_to_plot = {'Msn2(WT|4E|A)','Msn2(WT|4E|WT)','Msn2(WT|4E|T)'};

% output_folder = fullfile(parent_folder,'light_sweep_experiments','paper_figures');
output_folder = model_solutions_folder;

%% Import model output
data_store = fileDatastore(fullfile(model_solutions_folder,'mCitrine_model_round_*.mat'),'ReadFcn',@load);
mCitrine_model = cell(size(data_store.Files,1),1);
idx = 1;

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

%% Define ideal pulses nuclear Msn2

close all
reporters_to_plot = {'glpT','TKL2','ALD3','CTT1','DCS2','HXK1','RTN2','SIP18','SIP18 D6','SIP18 A4','DDR2','HSP12'};

t_measured = unique(data_stats.time);
pulse_t = linspace(min(t_measured),max(t_measured),1000)';

Msn2_params_list = [
    1	0.00	0	0	0	1	1	1
    1	0.25	0	50	0	1	1	1
    1	0.50	0	50	0	1	1	1
    1	0.75	0	50	0	1	1	1
    1	1.00	0	50	0	1	1	1
    1	1.00	0	10	0	1	1	1
    1	1.00	0	20	0	1	1	1
    1	1.00	0	30	0	1	1	1
    1	1.00	0	40	0	1	1	1
    2	1.00	0	5	5	6	1	1
    
    ];

Msn2_params_list = array2table(Msn2_params_list);
Msn2_params_list.Properties.VariableNames = {'signal_type','A','t0','t1','t2','cycles','c1','c2'};

Msn2_pulses_ideal = cell(size(Msn2_params_list,1),12);
idx = 1;
for pulse_idx = 1:size(Msn2_params_list,1)

    Msn2_params = Msn2_params_list(pulse_idx,:);
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
    Msn2_pulses_ideal{idx,10} = pulse_t;
    Msn2_pulses_ideal{idx,11} = Msn2_CT(pulse_t,Msn2_params);
    Msn2_pulses_ideal{idx,12} = trapz(pulse_t,Msn2_CT(pulse_t,Msn2_params));

    idx = idx + 1;
end

Msn2_pulses_ideal = cell2table(Msn2_pulses_ideal,'VariableNames',...
    {'pulse_idx','signal_type','A','t0','t1','t2','cycles','c1','c2',...
    'pulse_t','pulse_y','mScarlet_AUC'});
Msn2_pulses_ideal.pulse_idx(:,1) = 1:size(Msn2_pulses_ideal,1);
pulse_idx_list = unique(Msn2_pulses_ideal.pulse_idx);

%% Calculate toy promoter response to ideal pulses of nuclear Msn2

initial_conditions = [1 0 0 0 0 0];
n = 1;
K_list = 2.^(0:8);
kinetic_param_scale_list = 4.^(-4:4);
kinetic_param_list = {'k1','d1','k2','d2','k3'};
param_list = promoter_params_LHS.Properties.VariableNames;

fraction_active = 1;
K_scale_list = [0.5, 1, 2];

toy_promoter_response = cell(numel(K_list)*...
    numel(kinetic_param_scale_list)*numel(kinetic_param_list)*size(Msn2_pulses_ideal,1),20);
idx = 1;
for K_scale_idx = 1:numel(K_scale_list)

    K_scale = K_scale_list(K_scale_idx);

    for K_idx = 1:numel(K_list)

        K = K_list(K_idx);

        for kinetic_param_idx = 1:numel(kinetic_param_list)

            param_to_scale = kinetic_param_list(kinetic_param_idx);
            param_to_scale_idx = ismember(param_list,param_to_scale);

            param_scale_matrix = nan(1,numel(param_list));
            param_scale_matrix(param_to_scale_idx) = 1;
            param_scale_matrix = kinetic_param_scale_list'*param_scale_matrix;
            param_scale_matrix(isnan(param_scale_matrix)) = 1;
            param_scale_matrix(:,5) = n;
            param_scale_matrix(:,4) = K;

            for pulse_idx_temp = 1:size(Msn2_pulses_ideal,1)

                t_temp = Msn2_pulses_ideal.pulse_t(pulse_idx_temp);
                t_temp = t_temp{:};
                Msn2_temp = Msn2_pulses_ideal.pulse_y(pulse_idx_temp);
                Msn2_temp = Msn2_temp{:};

                % Initialize variables
                P_unbound_model_temp = zeros(numel(t_measured),numel(kinetic_param_scale_list));
                P_bound_model_temp = zeros(numel(t_measured),numel(kinetic_param_scale_list));
                P_active_model_temp = zeros(numel(t_measured),numel(kinetic_param_scale_list));
                mRNA_model_temp = zeros(numel(t_measured),numel(kinetic_param_scale_list));
                mCitrine_model_temp = zeros(numel(t_measured),numel(kinetic_param_scale_list));

                %              parfor kinetic_param_scale_idx = 1:numel(kinetic_param_scale_list)
                parfor kinetic_param_scale_idx = 1:numel(kinetic_param_scale_list)

                    promoter_params_temp = param_scale_matrix(kinetic_param_scale_idx,:);

                    [~,y] = ode45(@(t,y) promoter_ODEs(t,y,t_temp,Msn2_temp,promoter_params_temp,K_scale,fraction_active),...
                        t_measured,initial_conditions);
                    P_unbound_model_temp(:,kinetic_param_scale_idx) = y(:,1);
                    P_bound_model_temp(:,kinetic_param_scale_idx) = y(:,2);
                    P_active_model_temp(:,kinetic_param_scale_idx) = y(:,3);
                    mRNA_model_temp(:,kinetic_param_scale_idx) = y(:,4);
                    mCitrine_model_temp(:,kinetic_param_scale_idx) = y(:,6);

                end

                for kinetic_param_scale_idx = 1:numel(kinetic_param_scale_list)

                    pulse_idx = Msn2_pulses_ideal.pulse_idx(pulse_idx_temp);
%                     A = Msn2_pulses_ideal.A(pulse_idx_temp);
%                     t0 = Msn2_pulses_ideal.t0(pulse_idx_temp);
%                     t1 = Msn2_pulses_ideal.t1(pulse_idx_temp);
%                     t2 = Msn2_pulses_ideal.t2(pulse_idx_temp);
%                     cycles = Msn2_pulses_ideal.cycles(pulse_idx_temp);

                    [mCitrine_max_temp, idx_mCitrine_max_temp] = max(mCitrine_model_temp(:,kinetic_param_scale_idx));

                    toy_promoter_response{idx,1} = pulse_idx;
                    toy_promoter_response{idx,2} = Msn2_pulses_ideal.A(pulse_idx);
                    toy_promoter_response{idx,3} = Msn2_pulses_ideal.t0(pulse_idx);
                    toy_promoter_response{idx,4} = Msn2_pulses_ideal.t1(pulse_idx);
                    toy_promoter_response{idx,5} = Msn2_pulses_ideal.t2(pulse_idx);
                    toy_promoter_response{idx,6} = Msn2_pulses_ideal.cycles(pulse_idx);

                    toy_promoter_response{idx,7} = param_to_scale;
                    toy_promoter_response{idx,8} = kinetic_param_scale_list(kinetic_param_scale_idx);
                    toy_promoter_response{idx,9} = param_scale_matrix(kinetic_param_scale_idx,:);
                    toy_promoter_response{idx,10} = K_scale;

                    toy_promoter_response{idx,11} = t_measured;
                    toy_promoter_response{idx,12} = interp1(t_temp,Msn2_temp,t_measured);
                    toy_promoter_response{idx,13} = trapz(t_temp,Msn2_temp);

                    toy_promoter_response{idx,14} = P_unbound_model_temp(:,kinetic_param_scale_idx);
                    toy_promoter_response{idx,15} = P_bound_model_temp(:,kinetic_param_scale_idx);
                    toy_promoter_response{idx,16} = P_active_model_temp(:,kinetic_param_scale_idx);
                    toy_promoter_response{idx,17} = mRNA_model_temp(:,kinetic_param_scale_idx);
                    toy_promoter_response{idx,18} = mCitrine_model_temp(:,kinetic_param_scale_idx);
                    toy_promoter_response{idx,19} = max(P_active_model_temp(:,kinetic_param_scale_idx));
                    toy_promoter_response{idx,20} = mCitrine_max_temp;        
     
                    idx = idx + 1;
                end
            end
        end
    end
end


toy_promoter_response = cell2table(toy_promoter_response,'VariableNames',...
    {'pulse_idx','A','t0','t1','t2','cycles',...
    'param_scaled','scale_factor','params','K_scale'...
    'time','mScarlet_localization','mScarlet_AUC',...
    'P_unbound','P_bound','P_active','mRNA','mCitrine','P_active_max','mCitrine_max'});
toy_promoter_response = splitvars(toy_promoter_response,'params','NewVariableNames',param_list);
toy_promoter_response.param_scaled = categorical(toy_promoter_response.param_scaled);

%% Calculate normalization factors

% Normalize mCitrine for given promoter at K = 1
toy_promoter_response_K_scale_1 = toy_promoter_response(toy_promoter_response.K==1 & ...
    toy_promoter_response.K_scale==1,{'pulse_idx','param_scaled','scale_factor','mCitrine_max'});
toy_promoter_response_K_scale_1.Properties.VariableNames('mCitrine_max') = {'mCitrine_max_K_scale_1'};
toy_promoter_response = outerjoin(toy_promoter_response,toy_promoter_response_K_scale_1,'Type', 'Left', 'MergeKeys', true);

toy_promoter_response_WT = toy_promoter_response(toy_promoter_response.K_scale==1,{'pulse_idx','param_scaled','scale_factor','K','mCitrine_max'});
toy_promoter_response_WT.Properties.VariableNames('mCitrine_max') = {'mCitrine_max_WT'};
toy_promoter_response = outerjoin(toy_promoter_response,toy_promoter_response_WT,'Type', 'Left', 'MergeKeys', true);

%% %%%%%%%%%%%%%%%%%%%%%%%% PLOT MAIN FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 3A: plot response vs baseline

close all

clc; clear g; figure('position',[100 100 800 260]);
g = gramm('x',log2(toy_promoter_response.scale_factor),'y',log2(toy_promoter_response.K),...
    'color',(toy_promoter_response.mCitrine_max./toy_promoter_response.mCitrine_max_K_scale_1),...
    'subset',toy_promoter_response.mCitrine_max>0  & toy_promoter_response.pulse_idx==5 & toy_promoter_response.K_scale==1);
g.facet_grid([],toy_promoter_response.param_scaled);
g.geom_point();
g.set_continuous_color();
g.set_continuous_color('CLim',[0 1]);
g.set_color_options('map','jet')
g.set_names('x','log_{2}(scale factor)','y','log_{2}(K)','column','','row','','color','');
g.set_order_options('column',kinetic_param_list);
g.set_point_options('base_size',19,'markers',{'s'});
g.set_text_options('Interpreter','tex');
g.set_title('mCitrine vs baseline')
g.axe_property('XTick',[-8 -4 0 4 8]);
% g.no_legend()
g.draw();
g.redraw(0.075)
% export_fig(fullfile(pwd,'heatmap_mCitrine_vs_baseline'),'-png','-m4');
% export_fig(fullfile(pwd,'heatmap_mCitrine_vs_baseline_colormap'),'-png','-m4');

%% Figure 3B: plot response w/ K scaling vs K
Msn2_to_plot = {'Msn2(WT|4E|A)','Msn2(WT|4E|WT)','Msn2(WT|4E|T)'};

subset = toy_promoter_response.k1==1 & toy_promoter_response.d1==1 & ...
    toy_promoter_response.k2==1 & toy_promoter_response.d2==1 & ...
    toy_promoter_response.k3==1 & toy_promoter_response.mCitrine_max>0 & ...
    toy_promoter_response.pulse_idx==5 ;

close all
clc; clear g; figure('position',[100 100 250 250]);
g = gramm('x',log2(toy_promoter_response.K),...
    'y',(toy_promoter_response.mCitrine_max./toy_promoter_response.mCitrine_max_WT),...
    'color',(toy_promoter_response.K_scale),...
    'subset',subset);
g.stat_summary('geom','bar','dodge',0,'width',0.75);
g.set_names('x','log_{2}(K)','y','rel. expression','column','','row','','color','');
g.set_order_options('column',kinetic_param_list);
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.set_text_options('Interpreter','tex');
g.set_title('mCitrine vs K')
g.axe_property('XTick',[0 2 4 6 8],'XLim',[-0.5 8.5],'YTick',[0 1 2]);
g.no_legend();
g.draw();
g.redraw(0.1)
% export_fig(fullfile(pwd,'heatmap_mCitrine_vs_K'),'-png','-m4');

