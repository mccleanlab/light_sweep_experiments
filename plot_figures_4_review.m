clearvars; clc; close all

%% Load promoter fits
% Load promoter fits calculated when using meausured nuclear Msn2 time
% courses as model input TF(t)
promoter_fits_measured = load('promoter_fits_measured.mat');
promoter_fits_measured = promoter_fits_measured.promoter_fits; 
promoter_fits_measured.label(:,1) = categorical("measured");

% Load promoter fits calculated when using ideal nuclear Msn2 time courses
% (w/o degradation) as model input TF(t)
promoter_fits_ideal = load('promoter_fits_ideal.mat');
promoter_fits_ideal = promoter_fits_ideal.promoter_fits;
promoter_fits_ideal.label(:,1) = categorical("ideal");

%% Normalize RSS error by minimim RSS error calculated for given reporter/Msn2 mutant
promoter_fits = [promoter_fits_measured; promoter_fits_ideal];
promoter_fits.RSS_norm_fc(:,1) = nan;

promoter_fits_RSS_min = grpstats(promoter_fits,{'strain','plasmid'},'min','DataVars','RSS_norm');
% promoter_fits_RSS_min = grpstats(promoter_fits,{'strain','plasmid'},'min','DataVars','RSS_norm');
promoter_fits = outerjoin(promoter_fits,promoter_fits_RSS_min,'Type', 'Left', 'MergeKeys', true);
promoter_fits.RSS_norm_fc = promoter_fits.RSS_norm./promoter_fits.min_RSS_norm;

%% Figure S3D: plot overall RSS error for measured vs ideal 
close all
reporters_to_plot = {'RTN2','TKL2','SIP18','ALD3','CTT1','SIP18 D6','DCS2','SIP18 A4','DDR2','HXK1','HSP12'};

clc; clear g; figure('position',[100 100 400 350]);
g = gramm('x',promoter_fits.label,'y',promoter_fits.RSS_norm_fc,...
    'subset',ismember(promoter_fits.reporter,reporters_to_plot) & ...
    promoter_fits.plasmid~='pMM0845' & promoter_fits.rank<=10);
% g.facet_wrap(promoter_fits.plasmid,'ncols',3,'scale','independent')
g.stat_summary('type','ci','geom',{'bar','black_errorbar'},'setylim',true);
g.set_names('x','','y',['RSS/RSS min' newline '(top 0.1% parameter sets)']);
g.axe_property('YLim',[1 1.2]);
g.set_color_options('map',[175 175 175]/255);
g.set_order_options('x',0)
g.draw();
export_fig(fullfile(pwd,'RSS_ideal_vs_measured'),'-png','-m4');