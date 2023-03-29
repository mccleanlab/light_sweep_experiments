close all; clearvars; clc;

% This script calculates RSS error between reporter expression measurements
% in response to Msn2* and predicted expression obtained by simulating
% promoter response to Msn2(A)* while scaling K and n. Do same again for
% Msn2(T)*

%% %%%%%%%%%%%%% LOAD PRECALCULATED VARIABLES (IF APPLICABLE) %%%%%%%%%%%%%

load('plot_figures_2_data.mat')
return

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IMPORT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set parameters
n_guesses_plot = 10;
initial_conditions = [1 0 0 0 0 0];
plasmids_to_calc_RSS = {'pMM0847','pMM1079'}; % Calculate RSS for measurements for these plasmids (Msn2* and Msn2(T)*)
plasmids_to_calc_RSS = categorical(plasmids_to_calc_RSS);

%% Import measurements, simulations, etc

%%%%%%%%%%%%%%%%%%%%%%% Select input/output folders %%%%%%%%%%%%%%%%%%%%%%%

parent_folder  = pwd;
if contains(parent_folder,'\')
    parent_folder_idx   = strfind(parent_folder,'\');
else
    parent_folder_idx   = strfind(parent_folder,'/');
end
parent_folder = parent_folder(1:parent_folder_idx(end)-1);

% %% Select input/output folders
% 
% % parent_folder = 'D:\Google Drive\light_sweep_shared';
% parent_folder = 'G:\My Drive\light_sweep_shared';
% 
% if ~isfolder(parent_folder)
%     parent_folder = 'D:\GoogleDriveUW\light_sweep_shared';
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Import measurements %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Import LHS promoter parameters
model_folder = fullfile(parent_folder,'promoter_model');
% model_solutions_folder = fullfile(model_folder,'output_1_Hill_1_Msn2_100k');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%% Import model output %%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%% %%%%%%%%%%% FIT MEAUREMENTS TO ODEs AND CALCULATE RSS ERROR %%%%%%%%%%%%
clc

% Set parameters
n_guesses = size(mCitrine_model,1);
t_measured = unique(data_stats.time);
condition_list = unique(data_stats.condition);
param_list = promoter_params_LHS.Properties.VariableNames;

% Organize Msn2 vs time
Msn2_measured_all = zeros(numel(t_measured),numel(condition_list));
for condition = 1:numel(condition_list)
    Msn2_measured_all(:,condition) = data_Msn2.mScarlet_localization(data_Msn2.condition==condition);
end

strain_list = unique(data_stats.strain,'stable');
[strain_list, strain_order] = sort(string(strain_list));
strain_list = categorical(strain_list);

reporter_list = unique(data_stats.reporter,'stable');
reporter_list = reporter_list(strain_order);
reporter_list = categorical(reporter_list);

plasmid_list = unique(data_stats.plasmid);
plasmid_list = categorical(sort(string(plasmid_list)));

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
        promoter_fits_temp.RSS_norm_sp(1:n_guesses,1) = RSS_norm_sp(1:n_guesses,:);
        promoter_fits_temp.RSS_norm_sp_fc(1:n_guesses,1) = RSS_norm_sp(1:n_guesses,:)/min(RSS_norm_sp);
        
        promoter_fits{idx,1} = promoter_fits_temp;
        idx = idx + 1;
        
    end
end

promoter_fits = vertcat(promoter_fits{:});

%% Calculate promoter response for measured nuclear Msn2(A)* while scaling K and n


fraction_active = 1;
K_scale_list = [0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 4, 6];
n_scale_list = [0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5];

reporters_to_plot = reporter_list;

promoter_response_measured_K_scale = cell(numel(reporters_to_plot)*n_guesses_plot*numel(K_scale_list)*numel(n_scale_list)*numel(condition_list),20);
idx = 1;
for ii = 1:numel(reporters_to_plot)
    reporter = categorical(reporters_to_plot(ii));
    strain = unique(data_stats.strain(data_stats.reporter==reporter));
    disp(strain)
    
    plasmid = categorical("pMM0846");
    Msn2 = unique(data_stats.Msn2(data_stats.plasmid==plasmid));
    Msn2_tex = unique(data_stats.Msn2_tex(data_stats.plasmid==plasmid));
    
    % Get top params for strain/plasmid where # conditions w/in bounds >= (max # conditions w/in bounds - 1)
    subset = promoter_fits.strain==strain & promoter_fits.plasmid==plasmid;
    promoter_params_strain = promoter_fits{subset,param_list};
    promoter_params_strain = promoter_params_strain(1:n_guesses_plot,:);
    
    for guess = 1:n_guesses_plot
        
        promoter_params_temp = promoter_params_strain(guess,:);
        
        for K_scale_idx = 1:numel(K_scale_list)
            
            K_scale = K_scale_list(K_scale_idx);
            
            for n_scale_idx = 1:numel(n_scale_list)
                
                n_scale = n_scale_list(n_scale_idx);
                
                P_unbound_model_temp = zeros(numel(t_measured),numel(condition_list));
                P_bound_model_temp = zeros(numel(t_measured),numel(condition_list));
                P_active_model_temp = zeros(numel(t_measured),numel(condition_list));
                mRNA_model_temp = zeros(numel(t_measured),numel(condition_list));
                mCitrine_model_temp = zeros(numel(t_measured),numel(condition_list));
                
                parfor condition = 1:numel(condition_list)
                    Msn2_measured_temp = Msn2_measured_all(:,condition);
                    
                    [~,y] = ode45(@(t,y) promoter_ODEs_scale_K_and_n(t,y,t_measured,Msn2_measured_temp,promoter_params_temp,K_scale,n_scale,fraction_active),...
                        t_measured,initial_conditions);
                    P_unbound_model_temp(:,condition) = y(:,1);
                    P_bound_model_temp(:,condition) = y(:,2);
                    P_active_model_temp(:,condition) = y(:,3);
                    mRNA_model_temp(:,condition) = y(:,4);
                    mCitrine_model_temp(:,condition) = y(:,6);
                end
                
                for condition = 1:numel(condition_list)
                    promoter_response_measured_K_scale{idx,1} = strain;
                    promoter_response_measured_K_scale{idx,2} = reporter;
                    promoter_response_measured_K_scale{idx,3} = plasmid;
                    promoter_response_measured_K_scale{idx,4} = Msn2;
                    promoter_response_measured_K_scale{idx,5} = Msn2_tex;
                    
                    promoter_response_measured_K_scale{idx,6} = condition;
                    promoter_response_measured_K_scale{idx,7} = t_measured;
                    promoter_response_measured_K_scale{idx,8} = Msn2_measured_all(:,condition);
                    promoter_response_measured_K_scale{idx,9} = trapz(t_measured,Msn2_measured_all(:,condition));
                    
                    promoter_response_measured_K_scale{idx,10} = guess;
                    promoter_response_measured_K_scale{idx,11} = promoter_params_temp;
                    promoter_response_measured_K_scale{idx,12} = K_scale;
                    promoter_response_measured_K_scale{idx,13} = fraction_active;
                    
                    promoter_response_measured_K_scale{idx,14} = P_unbound_model_temp(:,condition);
                    promoter_response_measured_K_scale{idx,15} = P_bound_model_temp(:,condition);
                    promoter_response_measured_K_scale{idx,16} = P_active_model_temp(:,condition);
                    promoter_response_measured_K_scale{idx,17} = mRNA_model_temp(:,condition);
                    promoter_response_measured_K_scale{idx,18} = mCitrine_model_temp(:,condition);
                    promoter_response_measured_K_scale{idx,19} = max(mCitrine_model_temp(:,condition));
                    
                    promoter_response_measured_K_scale{idx,20} = n_scale;
                    
                    idx = idx + 1;
                end
            end
        end
    end
end

promoter_response_measured_K_scale = cell2table(promoter_response_measured_K_scale,'VariableNames',...
    {'strain','reporter','plasmid','Msn2','Msn2_tex',...
    'condition','time','mScarlet_nuclear','mScarlet_AUC',...
    'guess','params','K_scale','fraction_active',...
    'P_unbound','P_bound','P_active','mRNA','mCitrine','mCitrine_max',...
    'n_scale'});

%% Calculate RSS vs K_scale and n_scale



promoter_K_scale_RSS = cell(numel(reporters_to_plot)*numel(K_scale_list)*numel(n_scale_list)*numel(plasmids_to_calc_RSS),8);

idx = 1;
for ii = 1:numel(reporters_to_plot)
    reporter = categorical(reporters_to_plot(ii));
    strain = unique(data_stats.strain(data_stats.reporter==reporter));
    
    for K_scale_idx = 1:numel(K_scale_list)
        K_scale = K_scale_list(K_scale_idx);
        
        for n_scale_idx = 1:numel(n_scale_list)
            n_scale = n_scale_list(n_scale_idx);
            
            subset = promoter_response_measured_K_scale.strain==strain & ...
                promoter_response_measured_K_scale.plasmid=='pMM0846' & ...
                promoter_response_measured_K_scale.condition<=9 & ...
                promoter_response_measured_K_scale.K_scale==K_scale & ...
                promoter_response_measured_K_scale.n_scale==n_scale;
            vars_subset = {'strain','plasmid','condition','guess','mCitrine_max'};
            mCitrine_max_model_temp = promoter_response_measured_K_scale(subset,vars_subset);
            mCitrine_max_model_temp = unstack(mCitrine_max_model_temp,'mCitrine_max','guess');
            mCitrine_max_model_temp = mCitrine_max_model_temp{:,4:end}';
            
            for plasmid_idx = 1:numel(plasmids_to_calc_RSS)
                
                plasmid = plasmids_to_calc_RSS(plasmid_idx);
                Msn2 = unique(data_stats.Msn2(data_stats.plasmid==plasmid));
                Msn2_tex = unique(data_stats.Msn2_tex(data_stats.plasmid==plasmid));
                
                subset = mCitrine_stats.strain==strain & ...
                    mCitrine_stats.plasmid==plasmid & ...
                    mCitrine_stats.condition<=9;
                vars_subset = {'strain','plasmid','condition','replicate','mCitrine_max'};
                mCitrine_max_measured_temp = mCitrine_stats(subset,vars_subset);
                mCitrine_max_measured_temp = unstack(mCitrine_max_measured_temp,'mCitrine_max','replicate');
                mCitrine_max_measured_temp = mCitrine_max_measured_temp{:,4:end};
                mCitrine_max_measured_temp = reshape(mCitrine_max_measured_temp,1,size(mCitrine_max_measured_temp,1),size(mCitrine_max_measured_temp,2));
                
                RSS_temp = (mCitrine_max_measured_temp - mCitrine_max_model_temp).^2;
                RSS_temp = nansum(RSS_temp,'all');
                
                promoter_K_scale_RSS{idx,1} = strain;
                promoter_K_scale_RSS{idx,2} = reporter;
                promoter_K_scale_RSS{idx,3} = plasmid;
                promoter_K_scale_RSS{idx,4} = Msn2;
                promoter_K_scale_RSS{idx,5} = Msn2_tex;
                
                promoter_K_scale_RSS{idx,6} = K_scale;
                promoter_K_scale_RSS{idx,7} = RSS_temp;
                
                promoter_K_scale_RSS{idx,8} = n_scale;
                
                idx = idx + 1;
            end
        end
    end
end

promoter_K_scale_RSS = cell2table(promoter_K_scale_RSS,'VariableNames',...
    {'strain','reporter','plasmid','Msn2','Msn2_tex','K_scale','RSS','n_scale'});

promoter_K_scale_RSS_min = grpstats(promoter_K_scale_RSS,{'strain'},'min','DataVars','RSS');
promoter_K_scale_RSS_min = clean_grpstats(promoter_K_scale_RSS_min);
promoter_K_scale_RSS_min.Properties.VariableNames(end) = {'RSS_min'};
promoter_K_scale_RSS = outerjoin(promoter_K_scale_RSS,promoter_K_scale_RSS_min, 'Type', 'Left', 'MergeKeys', true);

%% %%%%%%%%%%%%%%%% SAVE VARIABLES TO AVOID RECALCULATING %%%%%%%%%%%%%%%%%

save('plot_figures_2_data',...
    'data_stats','data_Msn2','mCitrine_stats','promoter_fits',...
    'promoter_response_measured_K_scale','promoter_K_scale_RSS','promoter_K_scale_RSS_min');
return

%% %%%%%%%%%%%%%%%%%%%%%%%% PLOT MAIN FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure S4B: plot RSS/min(RSS) vs K_scale and n_scale 
close all
% Plot RSS/RSS_min vs K_scale (all reporters)
Msn2_to_plot = {'Msn2(WT|4E|WT)'};
clc; clear g; figure('position',[100 100 250 200]);
g = gramm('x',promoter_K_scale_RSS.K_scale,...
    'y',log(promoter_K_scale_RSS.RSS./promoter_K_scale_RSS.RSS_min),...
    'color',promoter_K_scale_RSS.n_scale,...
    'subset',~ismember(promoter_K_scale_RSS.reporter,'glpT') & ismember(promoter_K_scale_RSS.Msn2,Msn2_to_plot) & ...
    ismember(promoter_K_scale_RSS.n_scale,[0.5 1 1.5 1.75 2]));
% g.facet_wrap(promoter_K_scale_RSS.n_scale,'ncols',6);
g.stat_summary('type','std','geom','line','setylim',true);
% g.stat_summary('type','std','setylim',true);
g.set_names('x','','y','log(RSS/RSS_{min})','color','');
g.set_text_options('interpreter','tex');
g.axe_property('XLim',[0.5 4],'YLim',[0 12]);
g.no_legend()
g.draw();
g.redraw(0.05);
% export_fig(fullfile(pwd,'Msn2_min_RSS_scale_K'),'-png','-m4');

Msn2_to_plot = {'Msn2(WT|4E|T)'};
clc; clear g; figure('position',[100 100 250 200]);
g = gramm('x',promoter_K_scale_RSS.K_scale,...
    'y',log(promoter_K_scale_RSS.RSS./promoter_K_scale_RSS.RSS_min),...
    'color',promoter_K_scale_RSS.n_scale,...
    'subset',~ismember(promoter_K_scale_RSS.reporter,'glpT') & ismember(promoter_K_scale_RSS.Msn2,Msn2_to_plot) & ...
    ismember(promoter_K_scale_RSS.n_scale,[0.5 1 1.5 1.75 2]));
% g.facet_wrap(promoter_K_scale_RSS.n_scale,'ncols',6);
g.stat_summary('type','std','geom','line','setylim',true);
% g.stat_summary('type','std','setylim',true);
g.set_names('x','','y','log(RSS/RSS_{min})','color','');
g.set_text_options('interpreter','tex');
g.axe_property('XLim',[0.5 4],'YLim',[0 12]);
g.no_legend()
g.draw();
g.redraw(0.05);
% export_fig(fullfile(pwd,'Msn2T_min_RSS_scale_K'),'-png','-m4');

% Calculte K_scale value that minimizes RSS/RSS_min
% [~,RSS_min_idx] = min(g.results.stat_summary.y);
% K_scale_min = g.results.stat_summary.x(RSS_min_idx) ;

%% Figure S4C: plot simulated vs measured mCitrine for amplitude conditions

% % Plot RSS/RSS_min vs K_scale (all reporters)
% Msn2_to_plot = {'Msn2(WT|4E|WT)'};
% % Msn2_to_plot = {'Msn2(WT|4E|T)'};
% clc; clear g; figure('position',[100 100 400 350]);
% g = gramm('x',promoter_K_scale_RSS.K_scale,...
%     'y',log(promoter_K_scale_RSS.RSS./promoter_K_scale_RSS.RSS_min),...
%     'color',promoter_K_scale_RSS.n_scale,...
%     'subset',~ismember(promoter_K_scale_RSS.reporter,'glpT') & ismember(promoter_K_scale_RSS.Msn2,Msn2_to_plot) & ...
%     ismember(promoter_K_scale_RSS.n_scale,[0.5 1 1.5 1.75 2]));
% % g.facet_wrap(promoter_K_scale_RSS.n_scale,'ncols',6);
% g.stat_summary('type','std','geom','line','setylim',true);
% % g.stat_summary('type','std','setylim',true);
% g.set_names('x','','y','log(RSS/RSS_{min})','color','');
% % g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
% % g.set_order_options('color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
% g.set_text_options('interpreter','tex','base_size',12);
% g.axe_property('XLim',[0.5 4],'YLim',[0 12]);
% g.no_legend()
% % g.set_title('RSS/RSS_{min} vs K_{scale} for top 10 fits');
% g.draw();

% Calculate K_scale value that minimizes RSS/RSS_min
% [~,RSS_min_idx] = min(g.results.stat_summary.y);
% K_scale_min = g.results.stat_summary.x(RSS_min_idx) ;

%Label selected parameter sets
promoter_response_measured_K_scale.label(:,1) = categorical("NA");
promoter_response_measured_K_scale.label(promoter_response_measured_K_scale.K_scale==1 & ...
    promoter_response_measured_K_scale.n_scale==1)=categorical("Msn2(A)*: w = 1, m = 1");
promoter_response_measured_K_scale.label(promoter_response_measured_K_scale.K_scale==1.75 & ...
    promoter_response_measured_K_scale.n_scale==1)=categorical("Msn2*: w = 1.75, m = 1");
promoter_response_measured_K_scale.label(promoter_response_measured_K_scale.K_scale==1.75 & ...
    promoter_response_measured_K_scale.n_scale==1.75)=categorical("Msn2(T)*: w = 1.75, m = 1.75");

% Plot Figure 5A right panel: max mCitrine vs Msn2 AUC for K scaling of Msn2(WT|4E|A) (amplitude conditions)
close all
amplitude_conditions = 1:2:9;
duration_conditions = [1,2:2:8,9];
Msn2_to_plot ={'Msn2(WT|4E|A)','Msn2(WT|4E|WT)','Msn2(WT|4E|T)'};
% reporters_to_plot = {'glpT','TKL2','ALD3','SIP18 D6','DCS2','HXK1','RTN2','SIP18','CTT1','SIP18 A4','DDR2','HSP12'};
reporters_to_plot = {'glpT','RTN2','TKL2','SIP18','ALD3','CTT1','SIP18 D6','DCS2','SIP18 A4','DDR2','HXK1','HSP12'};
% reporters_to_plot = {'RTN2','CTT1','HSP12'};

% Plot measurements for amplitude conditions
% clc; clear g; figure('position',[100 100 800 500]);
clc; clear g; figure('position',[100 100 900 350]);
g = gramm('x',mCitrine_stats.mScarlet_AUC,'y',mCitrine_stats.mCitrine_max - mCitrine_stats.mCitrine_basal,...
    'color',mCitrine_stats.Msn2_tex,...
    'subset',ismember(mCitrine_stats.reporter,reporters_to_plot) & ismember(mCitrine_stats.Msn2,Msn2_to_plot) &...
    ismember(mCitrine_stats.condition,amplitude_conditions));
g.facet_wrap(mCitrine_stats.reporter,'ncols',6,'scale','independent');
g.stat_summary('type','std','geom',{'point','errorbar'});
% g.set_names('x','Msn2 AUC','y','max mCitrine','column','','row','','color','','linestyle','');
g.set_names('x','Msn2 AUC','y','','column','','row','','color','','linestyle','');
g.set_text_options('interpreter','tex');
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.set_order_options('column',reporters_to_plot,'color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.no_legend();
g.draw();

% Plot simulated for amplitude conditions
g.update('x',promoter_response_measured_K_scale.mScarlet_AUC,'y',promoter_response_measured_K_scale.mCitrine_max,...
    'color',promoter_response_measured_K_scale.label,...
    'subset',ismember(promoter_response_measured_K_scale.reporter,reporters_to_plot) & ...
    ismember(promoter_response_measured_K_scale.Msn2,Msn2_to_plot) & ...
    ismember(promoter_response_measured_K_scale.condition,amplitude_conditions) & promoter_response_measured_K_scale.label~='NA');
g.facet_wrap(promoter_response_measured_K_scale.reporter,'ncols',6,'scale','independent');
g.stat_summary('type','std');
g.set_order_options('column',reporters_to_plot,'color',0,'linestyle',-1);
g.set_line_options('base_size',1);
g.no_legend();
g.draw();
g.redraw(0.075);
% g.facet_axes_handles(1).YLim = [-1 20];
% export_fig(fullfile(pwd,'mCitrine_mutants_measured_vs_scaled_model'),'-png','-m4');



