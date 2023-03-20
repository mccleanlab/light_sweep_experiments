clearvars; clc; close all

%% Set folders and plot settings

% Set parent folder
% parent_folder = 'D:\Google Drive\light_sweep_shared\light_sweep_experiments';
% parent_folder = 'G:\My Drive\light_sweep_shared\light_sweep_experiments';
parent_folder = 'D:\GoogleDriveUW\light_sweep_shared\light_sweep_experiments';

% Set folder for data import
import_folder = fullfile(parent_folder,'collected_measurements');
% general_info_folder = 'D:\Google Drive\light_sweep_shared';
% general_info_folder = 'G:\My Drive\light_sweep_shared';
general_info_folder = 'D:\GoogleDriveUW\light_sweep_shared';

% Set output folder for graphs
output_folder = fullfile(parent_folder,'graphs_temp');

% Get hard-coded y limits for plots
opts = detectImportOptions(fullfile(general_info_folder,'plot_settings.xlsx'),'Sheet','mCitrine_y_limits');
opts.VariableTypes = {'categorical','double','double','double','double'};
mCitrine_y_limits = readtable(fullfile(general_info_folder,'plot_settings.xlsx'),opts);

% Get light on/off info
% general_info_folder = 'D:\Google Drive\light_sweep_shared';
opts = detectImportOptions(fullfile(general_info_folder,'plot_settings.xlsx'),'Sheet','t_light');
opts.VariableTypes = {'double','double','double'};
t_light = readtable(fullfile(general_info_folder,'plot_settings.xlsx'),opts);

% Get plot colors and Msn2_tex labels
% general_info_folder = 'D:\Google Drive\light_sweep_shared';
opts = detectImportOptions(fullfile(general_info_folder,'plot_settings.xlsx'),'Sheet','Msn2_tex');
opts.VariableTypes = {'categorical','categorical','double','double','double'};
plot_colors = readtable(fullfile(general_info_folder,'plot_settings.xlsx'),opts);
plot_colors.RGB = [plot_colors.R, plot_colors.G, plot_colors.B];
plot_colors = plot_colors(:,{'Msn2','Msn2_tex','RGB'});
plot_colors.RGB(plot_colors.Msn2=='Msn2',:) = [70, 115, 190];
% Msn2_to_plot = {'Msn2(WT|4E|A)','Msn2(WT|4E|WT)','Msn2(WT|4E|T)','Msn2'};
Msn2_to_plot = {'Msn2(WT|4E|A)','Msn2(WT|4E|WT)','Msn2(WT|4E|T)'};

%% Import measurements from .mat files and update labels
variables_to_import = {
    'experiment','well','condition','batch_control_type','replicate',...
    'frame','time','drop_frame',...
    'strain','plasmid', 'Msn2','Msn2_tex', 'CLASP', 'reporter',...
    'mScarlet_BG','mScarlet_nuclear_median','mScarlet_cyto_median','mScarlet_cell_median',...
    'mCitrine_BG','mCitrine_cell_median',...
    'amplitudes', 'pulse_start_times','pulse_numbs', 'pulse_high_times', 'pulse_low_times',...
    };

% Import data
data = import_mat_table(import_folder,variables_to_import);

% Make reporter names less annoying
data.reporter(data.reporter=='pSIP18_A4') = 'pSIP18 A4';
data.reporter(data.reporter=='pSIP18_D6') = 'pSIP18 D6';
data.reporter(data.reporter~='glpT') = categorical(extractAfter(string(data.reporter(data.reporter~='glpT')),1));

% Rename Msn2 mutants for tex
data.Msn2_tex(data.Msn2_tex=='None') = '\color{black}None';
data.Msn2_tex(data.Msn2_tex=='Msn2') = '\color{black}Msn2';

% Make list of strains reporters
strain_list = unique(data.strain,'stable');
[strain_list, strain_order] = sort(string(strain_list));
strain_list = categorical(strain_list);

% Make list of reporters
reporter_list = unique(data.reporter,'stable');
reporter_list = reporter_list(strain_order);
reporter_list = categorical(reporter_list);

% Force pulses for conditions 1 & 16 to 0 (was 1 pulse with 0 amplitude)
data.pulse_numbs(data.condition==1) = 0;
data.pulse_numbs(data.condition==16 & data.batch_control_type=='dark') = 0;

%% Process measurements

t_start = -10;
t_end = 152.5;
t_step = 2.5;

% Delete measurements past t_end
data(data.time>t_end,:) = [];

% Bin data based on time
bin_edges = t_start:t_step:t_end;
data.bin = discretize(data.time,bin_edges);

% Replace measurements from dropped frame with nan
data.mScarlet_nuclear(data.drop_frame==true) = nan;
data.mScarlet_cyto(data.drop_frame==true) = nan;
data.mCitrine_cell(data.drop_frame==true) = nan;

% Import irradiance correction model
% irradiance_model_folder = 'D:\Google Drive\light_sweep_shared\misc_experiments\irradiance_correction';
% irradiance_model_folder = 'G:\My Drive\light_sweep_shared\misc_experiments\irradiance_correction';
irradiance_model_folder = 'D:\GoogleDriveUW\light_sweep_shared\misc_experiments\irradiance_correction';
load(fullfile(irradiance_model_folder,'irradiance_correction_model.mat'));
a = irradiance_correction_model.a;
b = irradiance_correction_model.b;
c = irradiance_correction_model.c;

% Calculate scaling factor to correct for Intensilight brightness variation
day_0 = '20210713';
data.day = categorical(string(regexp(string(data.experiment),'\d{8}','match')));
data.day = datenum(string(data.day),'yyyymmdd') - datenum(day_0,'yyyymmdd');
data.irradiance_scale = c./(a*exp(-b*data.day) + c);

% Import mCitrine photobleaching model
% photobleach_model_folder = 'D:\Google Drive\light_sweep_shared\misc_experiments\photobleach_correction';
% photobleach_model_folder = 'G:\My Drive\light_sweep_shared\misc_experiments\photobleach_correction';
photobleach_model_folder = 'D:\GoogleDriveUW\light_sweep_shared\misc_experiments\photobleach_correction';
load(fullfile(photobleach_model_folder,'mCitrine_PBC_model.mat'));
a = mCitrine_PBC_model.a;
b = mCitrine_PBC_model.b;
c = mCitrine_PBC_model.c;

% Calculate scaling factor to correct mCitrine photobleaching
data.mCitrine_PBC_scale = 1./(a*exp(-b*data.frame) + c);

% Calculate BG fluorescence values for mScarlet and mCitrine (noise floor)
fluorescence_BG = grpstats(data(data.drop_frame~=true,:),'experiment','nanmedian','DataVars',{'mScarlet_BG','mCitrine_BG'});
fluorescence_BG = clean_grpstats(fluorescence_BG);
fluorescence_BG.Properties.VariableNames = regexprep(fluorescence_BG.Properties.VariableNames,'nanmedian_','');

fluorescence_BG.group(:,1) = 1;
fluorescence_BG = grpstats(fluorescence_BG,'group','nanmedian','DataVars',{'mScarlet_BG','mCitrine_BG'});
fluorescence_BG.Properties.VariableNames = regexprep(fluorescence_BG.Properties.VariableNames,'nanmedian_','');
fluorescence_BG = fluorescence_BG(:,{'mScarlet_BG','mCitrine_BG'});

% Subtract noise floor from nuclear mScarlet measurements and scale
data.mScarlet_nuclear = data.mScarlet_nuclear_median - fluorescence_BG.mScarlet_BG;
data.mScarlet_nuclear(data.mScarlet_nuclear<0) = 0;
data.mScarlet_nuclear = data.mScarlet_nuclear.*data.irradiance_scale;

% Subtract noise floor from cytoplasmic mScarlet measurements and scale
data.mScarlet_cyto = data.mScarlet_cyto_median - fluorescence_BG.mScarlet_BG;
data.mScarlet_cyto(data.mScarlet_cyto<=0) = 1;
data.mScarlet_cyto = data.mScarlet_cyto.*data.irradiance_scale;

% Subtract noise floor from cytoplasmic mScarlet measurements and scale
data.mScarlet_cell = data.mScarlet_cell_median - fluorescence_BG.mScarlet_BG;
data.mScarlet_cell(data.mScarlet_cell<=0) = 1;
data.mScarlet_cell = data.mScarlet_cell.*data.irradiance_scale;

% Subtract noise floor from mCitrine measurements and scale
data.mCitrine_cell = data.mCitrine_cell_median - fluorescence_BG.mCitrine_BG;
data.mCitrine_cell(data.mCitrine_cell<1) = 1;
data.mCitrine_cell = data.mCitrine_cell.*data.irradiance_scale;
data.mCitrine_cell = data.mCitrine_cell.*data.mCitrine_PBC_scale;

% Calculate mScarlet localization
data.mScarlet_localization = data.mScarlet_nuclear./data.mScarlet_cyto;

% % Delete mScarlet measurement ourliers (min and max 0.5% of measurements)
% mScarlet_lower_bound = prctile(data.mScarlet_nuclear,0.5);
% mScarlet_upper_bound = prctile(data.mScarlet_nuclear,99.5);
% data(data.mScarlet_nuclear<mScarlet_lower_bound | data.mScarlet_nuclear>mScarlet_upper_bound,:) = [];
%
% % Delete mCitrine measurement ourliers (min and max 0.5% of measurements)
% mCitrine_lower_bound = prctile(data.mCitrine_cell,0.1);
% mCitrine_upper_bound = prctile(data.mCitrine_cell(ismember(data.reporter,{'pHXK1','pHSP12'})),99.9);
% data(data.mCitrine_cell<mCitrine_lower_bound | data.mCitrine_cell>mCitrine_upper_bound,:) = [];

%% Calculate summary statistics

% Calculate median and std of mScarlet and mCitrine signals per bin
grp_vars = {'experiment','well','bin'};
data_vars = {'time','mScarlet_nuclear','mScarlet_localization','mCitrine_cell'};
data_stats_simple = grpstats(data,grp_vars,{'nanmedian'},'DataVars',data_vars);
data_stats_simple = clean_grpstats(data_stats_simple,false);
data_stats_simple.Properties.VariableNames;
data_stats_simple.Properties.VariableNames('GroupCount')={'cell_count_frame'};
data_stats_simple.Properties.VariableNames = regexprep(data_stats_simple.Properties.VariableNames,'nanmedian_','');

grp_vars = {'experiment','well','bin'};
data_vars = {'mCitrine_cell'};
mCitrine_noise = grpstats(data,grp_vars,@(x) (nanstd(x)./nanmean(x)).^2, 'DataVars',data_vars);
mCitrine_noise = clean_grpstats(mCitrine_noise);
mCitrine_noise.Properties.VariableNames(end) = {'mCitrine_noise'};

% Calculate median time per bin
data_time = grpstats(data_stats_simple,{'bin'},'nanmedian','DataVars','time');
data_time = clean_grpstats(data_time);
data_time.Properties.VariableNames = regexprep(data_time.Properties.VariableNames,'nanmedian_','');
data_stats_simple.time = [];

% Interpolate measurements per replicate the painful way
experiment_list = unique(data.experiment);
num_wells = 32;
num_bins = numel(unique(data.bin));

% Create list of all possible bins per experiment per well
bin_list = cell(numel(experiment_list)*num_wells*num_bins,1);
idx = 1;
for ii = 1:numel(experiment_list)
    experiment = experiment_list(ii);
    well_list = unique(data.well(data.experiment==experiment),'stable');
    
    for jj = 1:numel(well_list)
        well = well_list(jj);
        data_temp = table('Size',[num_bins, 3],...
            'VariableTypes',{'categorical','categorical','double'},...
            'VariableNames',{'experiment','well','bin'});
        
        data_temp.experiment(:,1) = experiment;
        data_temp.well(:,1) = well;
        data_temp.bin(:,1) = (1:num_bins)';
        
        bin_list{idx} = data_temp;
        idx = idx + 1;
    end
end
bin_list = vertcat(bin_list{:});

% Merge measurements into bin list and interpolate
data_stats = outerjoin(bin_list,data_stats_simple, 'Type', 'Left', 'MergeKeys', true);
data_stats = grouptransform(data_stats,{'experiment','well'},'linearfill',{'mScarlet_nuclear','mScarlet_localization','mCitrine_cell'});
data_stats = grouptransform(data_stats,{'experiment','well'},@(x) movmean(x,5),'mCitrine_cell');

% Merge time info back into data_stats
data_stats = outerjoin(data_stats,data_time, 'Type', 'Left', 'MergeKeys', true);
data_stats = sortrows(data_stats,{'experiment','well','bin'},{'ascend','ascend','ascend'});

% Apply labels to data_stats
grp_vars = {'experiment','well','strain','reporter','plasmid','Msn2','Msn2_tex','CLASP',...
    'condition','batch_control_type','replicate',...
    'amplitudes', 'pulse_start_times','pulse_numbs', 'pulse_high_times', 'pulse_low_times'};
data_labels = unique(data(:,grp_vars),'rows');
data_stats = outerjoin(data_stats,data_labels, 'Type', 'Left', 'MergeKeys', true);

% Calculate basal fluorescence
grp_vars = {'experiment','plasmid'};
data_vars = {'mScarlet_nuclear','mScarlet_localization','mCitrine_cell'};
fluorescence_basal = grpstats(data_stats(data_stats.time<0,:),grp_vars,'nanmedian','DataVars',data_vars);
fluorescence_basal = clean_grpstats(fluorescence_basal);
fluorescence_basal.Properties.VariableNames(end-2:end) = {'mScarlet_nuclear_basal','mScarlet_localization_basal','mCitrine_cell_basal'};
data_stats = outerjoin(data_stats,fluorescence_basal, 'Type', 'Left', 'MergeKeys', true);
data = outerjoin(data,fluorescence_basal, 'Type', 'Left', 'MergeKeys', true);

% Subtract basal fluorescence from data_stats nuclear signal
data_stats.mScarlet_nuclear = data_stats.mScarlet_nuclear - data_stats.mScarlet_nuclear_basal;
data_stats.mScarlet_nuclear(data_stats.mScarlet_nuclear<0) = 0;

data_stats.mScarlet_localization_keep_basal = data_stats.mScarlet_localization;
data_stats.mScarlet_localization = data_stats.mScarlet_localization - data_stats.mScarlet_localization_basal;
data_stats.mScarlet_localization(data_stats.mScarlet_localization<0) = 0;

% Calculate intial cell count per well
grp_vars = {'experiment','well'};
data_vars = {'cell_count_frame'};
cell_count_initial = grpstats(data_stats(data_stats.time<0,:),grp_vars,'nanmedian','DataVars',data_vars);
cell_count_initial = clean_grpstats(cell_count_initial);
cell_count_initial.Properties.VariableNames(end) = {'cell_count_initial'};
data_stats = outerjoin(data_stats,cell_count_initial, 'Type', 'Left', 'MergeKeys', true);

data_stats = outerjoin(data_stats,mCitrine_noise, 'Type', 'Left', 'MergeKeys', true);

% Apply re-ordered condition labels
condition_reorder = table();
condition_reorder.condition(:,1) = (1:14)';
condition_reorder.condition_display(:,1) = [1,6,2,7,3,8,4,9,5,12,11,10,13,14]';
data_stats = outerjoin(data_stats,condition_reorder, 'Type', 'Left', 'MergeKeys', true);
data_stats = sortrows(data_stats,'condition_display');

% Calculate max mCitrine per experiment per condition
grp_vars = {'experiment','well','strain','reporter','plasmid','Msn2','Msn2_tex','CLASP',...
    'condition','batch_control_type','replicate','condition_display',...
    'amplitudes', 'pulse_start_times','pulse_numbs', 'pulse_high_times', 'pulse_low_times'};
data_vars = {'mCitrine_cell'};
mCitrine_stats = grpstats(data_stats,grp_vars,'nanmax','DataVars',data_vars);
mCitrine_stats = clean_grpstats(mCitrine_stats);
mCitrine_stats.Properties.VariableNames(end) = {'mCitrine_max'};
mCitrine_stats = outerjoin(mCitrine_stats,fluorescence_basal(:,{'experiment','plasmid','mCitrine_cell_basal'}), 'Type', 'Left', 'MergeKeys', true);
mCitrine_stats.Properties.VariableNames('mCitrine_cell_basal') = {'mCitrine_basal'};

%% Plot Msn2 vs time for all reporters
close all

% reporters_to_exclude = {'pHXK1','pHSP12','pDDR2'};
% reporters_to_exclude = {'HXK1','HSP12','DDR2'};
% Msn2_to_plot = {'Msn2(WT|4E|A)','Msn2(WT|4E|WT)','Msn2(WT|4E|T)','Msn2','Msn2(WT|4E|D)'};

% clc; clear g; figure('units','normalized','outerposition',[0 0 1 0.6]);
% g = gramm('x',data_stats.time,'y',data_stats.mScarlet_localization,...
%     'color',data_stats.Msn2_tex,...
%     'subset',data_stats.condition<=14 & ~ismember(data_stats.reporter,reporters_to_exclude) & ismember(data_stats.Msn2,Msn2_to_plot));
% g.facet_wrap(data_stats.condition,'ncols',7);
% g.stat_summary('type','bootci','setylim',true);
% g.set_names('x','time (min.)','y','Msn2 localization','column','','color','');
% g.set_order_options('color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
% g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
% g.set_text_options('base_size',7,'interpreter','tex');
% g.set_title('Msn2 localization (nuc/cyt)');
% g.draw();
% export_fig(fullfile(output_folder,'Msn2_localization_per_mutant'),'-png','-m4');
% 
% clc; clear g; figure('units','normalized','outerposition',[0 0 1 0.6]);
% g = gramm('x',data_stats.time,'y',data_stats.mScarlet_nuclear,...
%     'color',data_stats.Msn2_tex,...
%     'subset',data_stats.condition<=14 & ismember(data_stats.Msn2,Msn2_to_plot));
% g.facet_wrap(data_stats.condition,'ncols',7);
% g.stat_summary('type','bootci','setylim',true);
% g.set_names('x','time (min.)','y','nuclear Msn2','column','','color','');
% g.set_order_options('color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
% g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
% g.set_text_options('base_size',7,'interpreter','tex');
% g.set_title('nuclear Msn2');
% g.axe_property('YLim',[-5 25]);
% g.draw();
% export_fig(fullfile(output_folder,'Msn2_nuclear_per_mutant'),'-png','-m4');

%% Calculate composite Msn2 vs time (for promoter model inputs)

grp_vars = {'condition','condition_display','bin','time'};
subset = data_stats.condition<=14 & ismember(data_stats.Msn2,{'Msn2(WT|4E|A)','Msn2(WT|4E|WT)','Msn2(WT|4E|T)'});
data_Msn2 = grpstats(data_stats(subset,:),grp_vars,'nanmean','DataVars',{'mScarlet_nuclear','mScarlet_localization'});
data_Msn2 = clean_grpstats(data_Msn2);
data_Msn2.Properties.VariableNames(end-1:end) = {'mScarlet_nuclear','mScarlet_localization'};

data_Msn2 = outerjoin(data_Msn2,t_light, 'Type', 'Left', 'MergeKeys', true);
data_Msn2.mScarlet_nuclear(data_Msn2.time<data_Msn2.t_light_start | data_Msn2.time>data_Msn2.t_light_end) = 0;

% Normalize to max
% data_stats.mScarlet_nuclear = data_stats.mScarlet_nuclear./max(data_Msn2.mScarlet_nuclear);
% data_Msn2.mScarlet_nuclear = data_Msn2.mScarlet_nuclear./max(data_Msn2.mScarlet_nuclear);
data_stats.mScarlet_localization = data_stats.mScarlet_localization./max(data_Msn2.mScarlet_localization);
data_Msn2.mScarlet_localization = data_Msn2.mScarlet_localization./max(data_Msn2.mScarlet_localization);

% close all
clc; clear g; figure('units','normalized','outerposition',[0 0 0.3 1]);
g = gramm('x',data_Msn2.time,'y',data_Msn2.mScarlet_nuclear,'subset',data_Msn2.condition<=14);
g.facet_wrap(data_Msn2.condition,'ncols',2);
g.stat_summary();
g.set_names('x','time (min.)','y','nuclear Msn2','column','');
g.draw();

clc; clear g; figure('units','normalized','outerposition',[0 0 0.3 1]);
g = gramm('x',data_Msn2.time,'y',data_Msn2.mScarlet_localization,'subset',data_Msn2.condition<=14);
g.facet_wrap(data_Msn2.condition,'ncols',2);
g.stat_summary();
g.set_names('x','time (min.)','y','Msn2 localization','column','');
g.draw();

%% Calculate Msn2 AUC
condition_list = unique(data_Msn2.condition);
data_Msn2_AUC = table('Size',[14,3],'VariableTypes',{'double','double','double'},'VariableNames',{'condition','condition_display','mScarlet_AUC'});

for ii = 1:numel(condition_list)
    condition = condition_list(ii);
    
    time = data_Msn2.time(data_Msn2.condition==condition);
    mScarlet_AUC = data_Msn2.mScarlet_localization(data_Msn2.condition==condition);
    AUC = trapz(time,mScarlet_AUC);
    
    data_Msn2_AUC.condition(ii,1) = condition;
    data_Msn2_AUC.condition_display(ii,1) = condition_reorder.condition_display(condition_reorder.condition==condition);
    data_Msn2_AUC.mScarlet_AUC(ii,1) = AUC; 
    
end

data_Msn2_AUC = outerjoin(data_Msn2_AUC,condition_reorder, 'Type', 'Left', 'MergeKeys', true);
mCitrine_stats = outerjoin(mCitrine_stats,data_Msn2_AUC, 'Type', 'Left', 'MergeKeys', true);


%% Calculate mCitrine AUC

experiment_list = unique(data_stats.experiment);
for ii = 1:numel(experiment_list)
    experiment = experiment_list(ii);
    well_list = unique(data_stats.well(data_stats.experiment==experiment));
    for jj = 1:numel(well_list)
        well = well_list(jj);
        time = data_stats.time(data_stats.experiment==experiment & data_stats.well==well);
        mCitrine_cell = data_stats.mCitrine_cell(data_stats.experiment==experiment & data_stats.well==well);
        AUC = trapz(time,mCitrine_cell);
        mCitrine_stats.mCitrine_AUC(mCitrine_stats.experiment==experiment & mCitrine_stats.well==well) = AUC;
    end
end

%% Export data to avoid above calculations

save('data.mat','data','-v7.3');
save('data_stats.mat','data_stats');
save('data_Msn2.mat','data_Msn2');
save('data_Msn2_AUC.mat','data_Msn2_AUC');
save('mCitrine_stats.mat','mCitrine_stats');

%% Count timelapses (79 experiments w/ 2528 timelapses)
n_experiments = unique(data_stats(:,{'experiment'}),'rows');
n_experiments = size(n_experiments,1)

n_timelapse = unique(data_stats(:,{'experiment','well'}),'rows');
n_timelapse = size(n_timelapse,1)
