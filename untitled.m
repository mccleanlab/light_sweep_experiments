

%% Figure S2A
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
export_fig(fullfile(pwd,'light_programs'),'-png','-m4');

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
export_fig(fullfile(pwd,'Msn2_loc_all_mutants'),'-png','-m4');

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
    export_fig(fullfile(pwd,strcat(string(reporter),'_mCitrine_vs_time_with_fits')),'-png','-m4');
end

%%
