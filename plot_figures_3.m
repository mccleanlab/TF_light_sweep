close all; clearvars; clc;

% Calculate amp/dur thresh from Pactive for ideal pulses Msn2

%% Select input/output folders

% parent_folder = 'D:\Google Drive\light_sweep_shared';
parent_folder = 'G:\My Drive\light_sweep_shared\';

if ~isfolder(parent_folder)
    parent_folder = 'D:\GoogleDriveUW\light_sweep_shared';
end

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

% Import parameters for ideal Msn2 localization function
opts = detectImportOptions(fullfile(parent_folder,'plot_settings.xlsx'),'Sheet','Msn2_CT_params');
Msn2_params_list = readtable(fullfile(parent_folder,'plot_settings.xlsx'),opts);

% Get plot colors and Msn2_tex labels
general_info_folder = 'D:\Google Drive\light_sweep_shared';
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

%% %%%%%%%%%%%% Load pre-calculated variables (if applicable) %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('plot_figures_3_data.mat')
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Fit measurements to ODE solutions
clc

% Set parameters
n_guesses = size(mCitrine_model,1);
initial_conditions = [1 0 0 0 0 0];
w_Msn2 = 1;
fraction_active = 1;

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
reporter_list(reporter_list=='glpT') = [];

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

%% Define ideal pulses nuclear Msn2

close all
reporters_to_plot = {'glpT','TKL2','ALD3','CTT1','DCS2','HXK1','RTN2','SIP18','SIP18 D6','SIP18 A4','DDR2','HSP12'};

pulse_end_list = 0:5:50;
pulse_amplitude_list = 0:0.10:1;
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
    Msn2_pulses_ideal{idx,3} = 1; % Msn2_params.A;
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

% Msn2_pulses_ideal = Msn2_pulses_ideal(11,:);
% Msn2_pulses_ideal = Msn2_pulses_ideal([1:2:11,19:2:25,34,36,37],:);
Msn2_pulses_ideal.A(Msn2_pulses_ideal.pulse_idx==1) = 0;
Msn2_pulses_ideal = Msn2_pulses_ideal(~ismember(Msn2_pulses_ideal.pulse_idx,[12,22]),:);

%% Calculate toy promoter response to ideal pulses of nuclear Msn2

K_list = 2.^(0:8);
% K_list = 2.^(0:2);
n_list = [0.5,1,1.5,2,4];
% n_list = [1,2];
kinetic_param_scale_list = 4.^(-4:4);
% kinetic_param_scale_list = 4.^(-1:1);
kinetic_param_list = {'k1','d1','k2','d2','k3'};

% n_param_scale_list = numel(kinetic_param_scale_list)*numel(kinetic_param_list);
toy_promoter_response = cell(numel(plasmid_list)*numel(n_list)*numel(K_list)*...
    numel(kinetic_param_scale_list)*numel(kinetic_param_list)*size(Msn2_pulses_ideal,1),27);
% toy_promoter_response = {};

idx = 1;
for plasmid_idx = 1:numel(plasmid_list)
    
    plasmid = plasmid_list(plasmid_idx);
    Msn2 = unique(data_stats.Msn2(data_stats.plasmid==plasmid));
    Msn2_tex = unique(data_stats.Msn2_tex(data_stats.plasmid==plasmid));
    
    if ismember(plasmid,'pMM0846')
        w_Msn2 = 1;
        m_Msn2 = 1;
    elseif ismember(plasmid,'pMM0847')
        w_Msn2 = 1.75;
        m_Msn2 = 1;
    elseif ismember(plasmid,'pMM1079')
        w_Msn2 = 1.75;
        m_Msn2 = 1.75;
    end
    
    disp(plasmid)

    for n_idx = 1:numel(n_list)
        
        n = n_list(n_idx);
        
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
                    
                    parfor kinetic_param_scale_idx = 1:numel(kinetic_param_scale_list)
                        
                        promoter_params_temp = param_scale_matrix(kinetic_param_scale_idx,:);
                        
                        [~,y] = ode45(@(t,y) promoter_ODEs_scale_K_and_n(t,y,t_temp,Msn2_temp,promoter_params_temp,w_Msn2,m_Msn2,fraction_active),...
                            t_measured,initial_conditions);
                        P_unbound_model_temp(:,kinetic_param_scale_idx) = y(:,1);
                        P_bound_model_temp(:,kinetic_param_scale_idx) = y(:,2);
                        P_active_model_temp(:,kinetic_param_scale_idx) = y(:,3);
                        mRNA_model_temp(:,kinetic_param_scale_idx) = y(:,4);
                        mCitrine_model_temp(:,kinetic_param_scale_idx) = y(:,6);
                        
                    end
                    
                    for kinetic_param_scale_idx = 1:numel(kinetic_param_scale_list)
                        
                        pulse_idx = Msn2_pulses_ideal.pulse_idx(pulse_idx_temp);
                        pulse_label = Msn2_pulses_ideal.pulse_label(pulse_idx_temp);
                        A = Msn2_pulses_ideal.A(pulse_idx_temp);
                        t0 = Msn2_pulses_ideal.t0(pulse_idx_temp);
                        t1 = Msn2_pulses_ideal.t1(pulse_idx_temp);
                        t2 = Msn2_pulses_ideal.t2(pulse_idx_temp);
                        cycles = Msn2_pulses_ideal.cycles(pulse_idx_temp);
                        
                        [mCitrine_max_temp, idx_mCitrine_max_temp] = max(mCitrine_model_temp(:,kinetic_param_scale_idx));
                        
                        toy_promoter_response{idx,1} = plasmid;
                        toy_promoter_response{idx,2} = Msn2;
                        toy_promoter_response{idx,3} = Msn2_tex;
                        
                        toy_promoter_response{idx,4} = pulse_idx;
                        toy_promoter_response{idx,5} = pulse_label;
                        toy_promoter_response{idx,6} = A;
                        toy_promoter_response{idx,7} = t0;
                        toy_promoter_response{idx,8} = t1;
                        toy_promoter_response{idx,9} = t2;
                        toy_promoter_response{idx,10} = cycles;
                        
                        toy_promoter_response{idx,11} = t_measured;
                        toy_promoter_response{idx,12} = interp1(t_temp,Msn2_temp,t_measured);
                        toy_promoter_response{idx,13} = trapz(t_temp,Msn2_temp);
                        
                        toy_promoter_response{idx,14} = param_to_scale;
                        toy_promoter_response{idx,15} = kinetic_param_scale_list(kinetic_param_scale_idx);
                        toy_promoter_response{idx,16} = param_scale_matrix(kinetic_param_scale_idx,:);
                        toy_promoter_response{idx,17} = w_Msn2;
                        toy_promoter_response{idx,18} = m_Msn2;
                        toy_promoter_response{idx,19} = fraction_active;
                        
                        toy_promoter_response{idx,20} = P_unbound_model_temp(:,kinetic_param_scale_idx);
                        toy_promoter_response{idx,21} = P_bound_model_temp(:,kinetic_param_scale_idx);
                        toy_promoter_response{idx,22} = P_active_model_temp(:,kinetic_param_scale_idx);
                        toy_promoter_response{idx,23} = mRNA_model_temp(:,kinetic_param_scale_idx);
                        toy_promoter_response{idx,24} = mCitrine_model_temp(:,kinetic_param_scale_idx);
                        toy_promoter_response{idx,25} = max(P_active_model_temp(:,kinetic_param_scale_idx));
                        toy_promoter_response{idx,26} = mCitrine_max_temp;
                        toy_promoter_response{idx,27} = t_measured(idx_mCitrine_max_temp);
                        
                        idx = idx + 1;
                    end
                end
            end
        end
    end
end
toy_promoter_response = cell2table(toy_promoter_response,'VariableNames',...
    {'plasmid','Msn2','Msn2_tex',...
    'pulse_idx','pulse_label','A','t0' 't1','t2','cycles',...
    'time','mScarlet_localization','mScarlet_AUC',...
    'param_scaled','scale_factor','params','w_Msn2','m_Msn2','fraction_active',...
    'P_unbound','P_bound','P_active','mRNA','mCitrine','P_active_max','mCitrine_max','t_mCitrine_max'});
toy_promoter_response = splitvars(toy_promoter_response,'params','NewVariableNames',param_list);
toy_promoter_response.param_scaled = categorical(toy_promoter_response.param_scaled);

%% Calculate additional parameters 
toy_promoter_response.activation_rate = (toy_promoter_response.k1.*toy_promoter_response.A.^toy_promoter_response.n)./...
    ((toy_promoter_response.w_Msn2.*toy_promoter_response.K.^toy_promoter_response.n) + toy_promoter_response.A.^toy_promoter_response.n);
toy_promoter_response.activation_ratio = toy_promoter_response.activation_rate./toy_promoter_response.d1;
toy_promoter_response.on_ratio = toy_promoter_response.k2./toy_promoter_response.d2;

% Msn2 localization model erroneously created small impulse of nuclear Msn2
% for zero amplitude pulse. Forcing predicted promoter activity and
% mCitrine expression to zero for these conditions
toy_promoter_response.P_active_max(toy_promoter_response.pulse_idx==1) = 0;
toy_promoter_response.mCitrine_max(toy_promoter_response.pulse_idx==1) = 0;
toy_promoter_response.mCitrine_max_vs_WT(toy_promoter_response.pulse_idx==1) = nan;

toy_promoter_response_param_idx = unique(toy_promoter_response(:,param_list),'rows');
toy_promoter_response_param_idx.param_idx(:,1) = 1:size(toy_promoter_response_param_idx,1);
toy_promoter_response = outerjoin(toy_promoter_response,toy_promoter_response_param_idx,'Type', 'Left', 'MergeKeys', true);

toy_promoter_response_full_light = toy_promoter_response(toy_promoter_response.pulse_idx==11,:);
toy_promoter_response_full_light = toy_promoter_response_full_light(:,{'Msn2','param_idx','mCitrine_max'});
toy_promoter_response_full_light.Properties.VariableNames(end) = {'mCitrine_max_full_light'};
toy_promoter_response = outerjoin(toy_promoter_response,toy_promoter_response_full_light,'Type', 'Left', 'MergeKeys', true);

toy_promoter_response_full_light_WT = toy_promoter_response(toy_promoter_response.Msn2=='Msn2(WT|4E|WT)' & toy_promoter_response.pulse_idx==11,:);
toy_promoter_response_full_light_WT = toy_promoter_response_full_light_WT(:,{'param_idx','mCitrine_max'});
toy_promoter_response_full_light_WT.Properties.VariableNames(end) = {'mCitrine_max_full_light_WT'};
toy_promoter_response = outerjoin(toy_promoter_response,toy_promoter_response_full_light_WT,'Type', 'Left', 'MergeKeys', true);

%% Calculate mCitrine vs WT for toy promoter
clc
toy_promoter_response.mCitrine_max_vs_WT(:,1) = nan;

% Plot fits over measurements
for param_idx = 1:numel(param_list)
    
    param_scaled = param_list(param_idx);
    disp(param_scaled)
    
    param_scale_list = unique(toy_promoter_response.scale_factor(ismember(toy_promoter_response.param_scaled,param_scaled)));
    pulse_idx_list_temp = unique(toy_promoter_response.pulse_idx);
    
    for n_idx = 1:numel(n_list)
        
        n = n_list(n_idx);
        
        for K_idx = 1:numel(K_list)
            
            K = K_list(K_idx);
            
            for kinetic_param_scale_idx = 1:numel(param_scale_list)
                
                scale_factor = param_scale_list(kinetic_param_scale_idx);
                
                for pulse_idx_temp = 1:numel(pulse_idx_list_temp)
                    
                    pulse_idx = pulse_idx_list_temp(pulse_idx_temp);
                    
                    subset = toy_promoter_response.pulse_idx==pulse_idx & ...
                        ismember(toy_promoter_response.param_scaled,param_scaled) & ...
                        toy_promoter_response.scale_factor==scale_factor & ...
                        toy_promoter_response.n==n & toy_promoter_response.K==K;
                    
                    toy_promoter_response_temp = toy_promoter_response(subset,:);
                    
                    mCitrine_max_WT = toy_promoter_response_temp.mCitrine_max(toy_promoter_response_temp.plasmid=='pMM0847');
                    mCitrine_max_vs_WT = toy_promoter_response_temp.mCitrine_max./mCitrine_max_WT;
                    
                    toy_promoter_response.mCitrine_max_vs_WT(subset,:) = mCitrine_max_vs_WT;
                    
                end
            end
        end
    end
end

%% Calculate slope ratios

continuous_conditions = [1,3,7,11];
pulsed_conditions = [1,25,24,23];
toy_promoter_thresholds = cell(numel(plasmid_list)*numel(n_list)*numel(K_list)*numel(kinetic_param_list)*numel(param_scale_list),1);

idx = 1;
for plasmid_idx = 1:numel(plasmid_list)
    
    plasmid = plasmid_list(plasmid_idx);
    
    disp(plasmid)
    
    for n_idx = 1:numel(n_list)
        
        n = n_list(n_idx);
        
        for K_idx = 1:numel(K_list)
            
            K = K_list(K_idx);
            
            for kinetic_param_idx = 1:numel(kinetic_param_list)
                
                param_scaled = kinetic_param_list(kinetic_param_idx);
                
                for scale_factor_idx = 1:numel(param_scale_list)
                    
                    scale_factor = param_scale_list(scale_factor_idx);
                    
                    %%%%%%%%%%%%%%%% Fit pulsed conditions %%%%%%%%%%%%%%%%
                    % Get predicted expression for continuous conditions
                    subset = (toy_promoter_response.plasmid==plasmid & ...
                        ismember(toy_promoter_response.pulse_idx,continuous_conditions) & ...
                        toy_promoter_response.n==n & ...
                        toy_promoter_response.K==K & ...
                        toy_promoter_response.param_scaled==param_scaled & ...
                        toy_promoter_response.scale_factor==scale_factor);
                    
                    toy_promoter_response_temp = toy_promoter_response(subset,:);
                    x_continuous = toy_promoter_response_temp.mScarlet_AUC;
                    y_continuous = toy_promoter_response_temp.mCitrine_max;
                    
                    % Fit predicted expression for continuous conditions
                    mdl_continuous = fitglm(x_continuous,y_continuous);
                    m_continuous = mdl_continuous.Coefficients.Estimate(2);
                    m_continuous_se = mdl_continuous.Coefficients.SE(2);
                    m_continuous_p = mdl_continuous.Coefficients.pValue(2);
                    b_continuous = mdl_continuous.Coefficients.Estimate(1);
                    
                    x_mdl_continuous = (linspace(min(x_continuous),max(x_continuous),100))';
                    y_mdl_continuous = m_continuous.*x_mdl_continuous + b_continuous;
                    
                    %%%%%%%%%%%%%%%% Fit pulsed conditions %%%%%%%%%%%%%%%%
                    % Get predicted expression for pulsed conditions
                    subset = (toy_promoter_response.plasmid==plasmid & ...
                        ismember(toy_promoter_response.pulse_idx,pulsed_conditions) & ...
                        toy_promoter_response.n==n & ...
                        toy_promoter_response.K==K & ...
                        toy_promoter_response.param_scaled==param_scaled & ...
                        toy_promoter_response.scale_factor==scale_factor);
                    
                    toy_promoter_response_temp = toy_promoter_response(subset,:);
                    x_pulsed = toy_promoter_response_temp.mScarlet_AUC;
                    y_pulsed = toy_promoter_response_temp.mCitrine_max;
                    
                    % Fit predicted expression for pulsed conditions
                    mdl_pulsed = fitglm(x_pulsed,y_pulsed);
                    m_pulsed = mdl_pulsed.Coefficients.Estimate(2);
                    m_pulsed_se = mdl_pulsed.Coefficients.SE(2);
                    m_pulsed_p = mdl_pulsed.Coefficients.pValue(2);
                    b_pulsed = mdl_pulsed.Coefficients.Estimate(1);
                    
                    x_mdl_pulsed = (linspace(min(x_pulsed),max(x_pulsed),100))';
                    y_mdl_pulsed = m_pulsed.*x_mdl_pulsed + b_pulsed;
                    
                    %%%%%%%%%%%% Plot to slopes to sanity check %%%%%%%%%%%
                    %                     close all; hold on
                    %                     plot(x_continuous,y_continuous,'ko')
                    %                     plot(x_mdl_continuous,y_mdl_continuous,'k-')
                    %                     plot(x_pulsed,y_pulsed,'k^')
                    %                     plot(x_mdl_pulsed,y_mdl_pulsed,'k:')
                    %
                    
                    %%%%%%%%%%%%%%%% Calculate slope ratio %%%%%%%%%%%%%%%%
                    slope_ratio = m_pulsed./m_continuous;
                    slope_ratio_se = (m_pulsed/m_continuous)* sqrt((m_pulsed_se/m_pulsed).^2 + (m_continuous_se/m_continuous).^2);
                                        
                    %%%%%%%%%%%%%%%%%%%%%%%%% Save %%%%%%%%%%%%%%%%%%%%%%%%
                    subset = (toy_promoter_response.plasmid==plasmid & ...
                        ismember(toy_promoter_response.pulse_idx,11) & ...
                        toy_promoter_response.n==n & ...
                        toy_promoter_response.K==K & ...
                        toy_promoter_response.param_scaled==param_scaled & ...
                        toy_promoter_response.scale_factor==scale_factor);
                    vars_to_keep = {'plasmid','Msn2','Msn2_tex',...
                        'param_scaled','scale_factor','param_idx',...
                        'k1','d1','k2','K','n','d2','k3',...
                        'w_Msn2','m_Msn2','fraction_active',...
                        'mCitrine_max'};
                    
                    toy_promoter_thresholds_temp = cell(1,5);
                    toy_promoter_thresholds_temp = toy_promoter_response(subset,vars_to_keep);
                    
                    toy_promoter_thresholds_temp.x_continuous = {x_continuous};
                    toy_promoter_thresholds_temp.y_continuous = {y_continuous};
                    toy_promoter_thresholds_temp.x_mdl_continuous = {x_mdl_continuous};
                    toy_promoter_thresholds_temp.y_mdl_continuous = {y_mdl_continuous};
                    
                    toy_promoter_thresholds_temp.x_pulsed = {x_pulsed};
                    toy_promoter_thresholds_temp.y_pulsed = {y_pulsed};
                    toy_promoter_thresholds_temp.x_mdl_pulsed = {x_mdl_pulsed};
                    toy_promoter_thresholds_temp.y_mdl_pulsed = {y_mdl_pulsed};
                    
                    toy_promoter_thresholds_temp.slope_continuous = m_continuous;
                    toy_promoter_thresholds_temp.slope_pulsed = m_pulsed;
                    toy_promoter_thresholds_temp.slope_ratio = slope_ratio;
                    toy_promoter_thresholds_temp.slope_ratio_se = slope_ratio_se;
                    
                    toy_promoter_thresholds{idx} = toy_promoter_thresholds_temp;    
                    
                    
                    idx = idx + 1;
                    
%                     if idx>10
%                         return
%                     end
                end
            end
        end
    end
end

toy_promoter_thresholds = vertcat(toy_promoter_thresholds{:});

toy_promoter_thresholds_WT = toy_promoter_thresholds(toy_promoter_thresholds.Msn2=='Msn2(WT|4E|WT)',{'param_idx','slope_ratio'});
toy_promoter_thresholds_WT.Properties.VariableNames(end) = {'slope_ratio_WT'};
toy_promoter_thresholds = outerjoin(toy_promoter_thresholds,toy_promoter_thresholds_WT,'Type', 'Left', 'MergeKeys', true);

%% Save variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('plot_figure_6_scale_n_data.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return

%% Figure 5B: heatmap of relative expression vs parameter space

% cmap = y./255;
cmap = viridis;
% cmap = circshift(cmap,125);
% cmap = cmap(20:235,:);
% cmap = (flipud(cmap))
% cmap = xx;
% cmap = xx/255;
% imshow(xx)
% cmap = rgb2gray(cmap); 
% return

clear g; clc; close all
clc; clear g; figure('position',[100 100 900 230]);
g = gramm('x',log2(toy_promoter_response.scale_factor),'y',log2(toy_promoter_response.K),...
    'color',(toy_promoter_response.mCitrine_max_vs_WT),...
    'subset',toy_promoter_response.Msn2=='Msn2(WT|4E|A)' & toy_promoter_response.n==1 & ...
    toy_promoter_response.mCitrine_max>0 & toy_promoter_response.pulse_idx==11);
g.facet_grid([],toy_promoter_response.param_scaled);
g.geom_point();
g.set_continuous_color('CLim',[0 2]);
% g.set_color_options('map','jet')
g.set_names('x','log_{2}(scale factor)','y','log_{2}(K)','column','','row','','color','');
g.set_order_options('column',kinetic_param_list);
g.set_point_options('base_size',15.5,'markers',{'s'});
g.set_text_options('Interpreter','tex');
g.axe_property('XTick',[-8 -4 0 4 8])
g.draw();
g.redraw(0.075)
% colormap(cmap)
export_fig(fullfile(pwd,'heatmap_params_A_vs_WT'),'-png','-m4');

% clear g; clc; close all
clc; clear g; figure('position',[100 100 900 230]);
g = gramm('x',log2(toy_promoter_response.scale_factor),'y',log2(toy_promoter_response.K),...
    'color',(toy_promoter_response.mCitrine_max_vs_WT),...
    'subset',toy_promoter_response.Msn2=='Msn2(WT|4E|T)' & toy_promoter_response.n==1 & ...
    toy_promoter_response.mCitrine_max>0 & toy_promoter_response.pulse_idx==11);
g.facet_grid([],toy_promoter_response.param_scaled);
g.geom_point();
% g.set_continuous_color('CLim',[0 1]);
g.set_continuous_color('CLim',[0 2]);
g.set_color_options('map','jet')
g.set_names('x','log_{2}(scale factor)','y','log_{2}(K)','column','','row','','color','');
g.set_order_options('column',kinetic_param_list);
g.set_point_options('base_size',15.5,'markers',{'s'});
g.set_text_options('Interpreter','tex');
g.axe_property('XTick',[-8 -4 0 4 8])
g.draw();
g.redraw(0.075)
% colormap(cmap)
export_fig(fullfile(pwd,'heatmap_params_T_vs_WT'),'-png','-m4');

%% Figure 5C: mCitrine vs pulse amplitude

close all
Msn2_to_plot = {'Msn2(WT|4E|A)','Msn2(WT|4E|WT)','Msn2(WT|4E|T)'};
param_order = {'k1','d1','K','n','k2','d2','k3'};
interpulse_conditions = [7,24,26,27];
duration_conditions = 1:11;
amplitude_conditions = [1,11,13:21];

toy_promoter_response.label(:,1) = categorical("NA");
toy_promoter_response.label(toy_promoter_response.K==2 & toy_promoter_response.n==1) = categorical("high");
toy_promoter_response.label(toy_promoter_response.K==4 & toy_promoter_response.n==2) = categorical("low");

promoter_label = 'high';
param_scaled = 'k1';
scale_factor = 16;

var_to_plot = 'A';
conditions_to_plot = amplitude_conditions;
x_label = 'amplitude (%)';

clc; clear g; figure('position',[100 100 250 250]);
g = gramm('x',100*(toy_promoter_response.(var_to_plot)),'y',(toy_promoter_response.P_active_max),...
    'color',toy_promoter_response.Msn2_tex,...
    'subset',toy_promoter_response.label==promoter_label & ...
    ismember(toy_promoter_response.param_scaled,param_scaled) & ...
    ismember(toy_promoter_response.scale_factor,scale_factor) & ...
    ismember(toy_promoter_response.pulse_idx,conditions_to_plot));
g.stat_summary();
g.set_names('x',x_label,'y','max P_{on}','column','','row','','color','')
g.set_order_options('color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.set_text_options('interpreter','tex','base_size',12);
g.no_legend();
% % g.set_line_options('styles',{':'});
g.draw();
export_fig(fullfile("fig_5_vary_amplitude_Pon"),'-png','-m4');

clc; clear g; figure('position',[100 100 250 250]);
g = gramm('x',100*(toy_promoter_response.(var_to_plot)),'y',(toy_promoter_response.mCitrine_max),...
    'color',toy_promoter_response.Msn2_tex,...
    'subset',toy_promoter_response.label==promoter_label & ...
    ismember(toy_promoter_response.param_scaled,param_scaled) & ...
    ismember(toy_promoter_response.scale_factor,scale_factor) & ...
    ismember(toy_promoter_response.pulse_idx,conditions_to_plot));
g.stat_summary();
g.set_names('x',x_label,'y','max mCitrine','column','','row','','color','')
g.set_order_options('color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.set_text_options('interpreter','tex','base_size',12);
g.no_legend();
g.draw();

g.update('subset',toy_promoter_response.label==promoter_label & ...
    ismember(toy_promoter_response.param_scaled,param_scaled) & ...
    ismember(toy_promoter_response.scale_factor,scale_factor) & ...
    ismember(toy_promoter_response.pulse_idx,conditions_to_plot) & ...
    ismember(toy_promoter_response.A,[0.2,0.5,0.8]));
g.geom_point();
g.no_legend();
g.draw();
g.redraw(0.025);
export_fig(fullfile("fig_5_vary_amplitude"),'-png','-m4');

clc; clear g; figure('position',[100 100 100 250]);
g = gramm('x',categorical(100*(toy_promoter_response.(var_to_plot))),'y',(toy_promoter_response.mCitrine_max_vs_WT),...
    'color',toy_promoter_response.Msn2_tex,...
    'subset',toy_promoter_response.label==promoter_label & ...
    ismember(toy_promoter_response.param_scaled,param_scaled) & ...
    ismember(toy_promoter_response.scale_factor,scale_factor) & ...
    ismember(toy_promoter_response.pulse_idx,conditions_to_plot) & ...
    ismember(toy_promoter_response.A,[0.2,0.5,0.8]));
g.set_names('x',x_label,'y','','column','','row','','color','')
g.geom_bar('width',0.5,'LineWidth',1)
g.set_order_options('x',0,'row',param_order,'color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.set_text_options('interpreter','tex','base_size',11);
g.no_legend();
g.draw()
set(findobj('color','g'),'Color',[200 200 200]/255);
g.redraw(0.025)
export_fig(fullfile("fig_5_vary_amplitude_relative"),'-png','-m4');

subset = mCitrine_stats.Msn2=='Msn2(WT|4E|WT)' & mCitrine_stats.condition<=14;
grp_vars = {'reporter','condition'}
mCitrine_stats_WT = grpstats(mCitrine_stats(subset,:),grp_vars,'nanmean','DataVars','mCitrine_max');
mCitrine_stats_WT = clean_grpstats(mCitrine_stats_WT);
mCitrine_stats_WT.Properties.VariableNames(end) = {'mCitrine_max_WT'}
mCitrine_stats = outerjoin(mCitrine_stats,mCitrine_stats_WT,'Type', 'Left', 'MergeKeys', true);

reporters_to_plot = 'DDR2';
duration_conditions = [1,2:2:8,9];
amplitude_conditions = 1:2:9;
clc; clear g; figure('position',[100 100 250 250]);
g = gramm('x',mCitrine_stats.amplitudes,'y',mCitrine_stats.mCitrine_max./mCitrine_stats.mCitrine_max_WT,...
    'color',mCitrine_stats.Msn2_tex,...
    'subset',ismember(mCitrine_stats.condition,amplitude_conditions) & ...
    ismember(mCitrine_stats.Msn2,Msn2_to_plot) & ismember(mCitrine_stats.reporter,reporters_to_plot));
g.facet_grid(mCitrine_stats.reporter,[],'row_labels',false,'scale','independent');
g.stat_summary('type','std','geom',{'point','errorbar'},'setylim',true,'width',4,'dodge',0.75);
g.set_names('x','Msn2 AUC','y','mCitrine','column','','row','','linestyle','','color','','marker','increasing duration)');
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.set_order_options('color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.set_text_options('base_size',10,'interpreter','tex','title_scaling',1);
g.set_point_options('base_size',6,'markers',{'o','^'});
% g.axe_property('XLim',[-5 55]);
g.no_legend();
g.draw()

%% Figure S12: plot heatmap of absolute expression vs parameter space
close all

clc; clear g; figure('position',[100 100 900 450]);
g = gramm('x',log2(toy_promoter_response.scale_factor),'y',log2(toy_promoter_response.K),...
    'color',log10(toy_promoter_response.mCitrine_max),...
    'subset',toy_promoter_response.Msn2=='Msn2(WT|4E|A)' & toy_promoter_response.pulse_idx==11);
g.facet_grid(toy_promoter_response.n,toy_promoter_response.param_scaled);
g.geom_point();
g.set_continuous_color('CLim',[0 5]);
g.set_names('x','log_{2}(scale factor)','y','log_{2}(K)','column','','row','','color','');
g.set_order_options('column',kinetic_param_list);
g.set_point_options('base_size',15,'markers',{'s'});
g.set_text_options('Interpreter','tex');
g.draw();
% g.redraw(0.075)

clc; clear g; figure('position',[100 100 900 450]);
g = gramm('x',log2(toy_promoter_response.scale_factor),'y',log2(toy_promoter_response.K),...
    'color',log10(toy_promoter_response.mCitrine_max),...
    'subset',toy_promoter_response.Msn2=='Msn2(WT|4E|WT)' & toy_promoter_response.pulse_idx==11);
g.facet_grid(toy_promoter_response.n,toy_promoter_response.param_scaled);
g.geom_point();
g.set_continuous_color('CLim',[0 5]);
g.set_names('x','log_{2}(scale factor)','y','log_{2}(K)','column','','row','','color','');
g.set_order_options('column',kinetic_param_list);
g.set_point_options('base_size',15,'markers',{'s'});
g.set_text_options('Interpreter','tex');
g.draw();
% g.redraw(0.075)

clc; clear g; figure('position',[100 100 900 450]);
g = gramm('x',log2(toy_promoter_response.scale_factor),'y',log2(toy_promoter_response.K),...
    'color',log10(toy_promoter_response.mCitrine_max),...
    'subset',toy_promoter_response.Msn2=='Msn2(WT|4E|T)' & toy_promoter_response.pulse_idx==11);
g.facet_grid(toy_promoter_response.n,toy_promoter_response.param_scaled);
g.geom_point();
g.set_continuous_color('CLim',[0 5]);
g.set_names('x','log_{2}(scale factor)','y','log_{2}(K)','column','','row','','color','');
g.set_order_options('column',kinetic_param_list);
g.set_point_options('base_size',15,'markers',{'s'});
g.set_text_options('Interpreter','tex');
g.draw();
% g.redraw(0.075)

%% Figure S13: plot expression vs amplitude, duration, interpulse vs parameter scaling
close all
Msn2_to_plot = {'Msn2(WT|4E|A)','Msn2(WT|4E|WT)','Msn2(WT|4E|T)'};
param_order = {'k1','d1','K','n','k2','d2','k3'};
interpulse_conditions = [7,24,26,27];
duration_conditions = 1:11;
amplitude_conditions = [1,13:21,11];

param_scaled = 'd2';
var_to_plot = 't2';
conditions_to_plot = interpulse_conditions;
x_label = 'interpulse duration (min)';
K_list = [1,4,64];
scale_factor_list = [1/16,1,16];
row_order = 0;

% clc; clear g; figure('position',[0 0 1000 500]);
% g = gramm('x',(toy_promoter_response.(var_to_plot)),'y',(toy_promoter_response.P_active_max),...
%     'color',toy_promoter_response.Msn2_tex,...
%     'subset',ismember(toy_promoter_response.pulse_idx,conditions_to_plot) & ...
%     toy_promoter_response.param_scaled==param_scaled & toy_promoter_response.n==1 & ...
%     ismember(toy_promoter_response.K,K_list) & ismember(toy_promoter_response.scale_factor,scale_factor_list));
% g.facet_grid(toy_promoter_response.scale_factor,toy_promoter_response.K,'scale','independent');
% % g.geom_line();
% g.stat_summary();
% g.set_names('x',x_label,'y','max P_{active}','column','','row','','color','')
% g.set_order_options('row',row_order,'color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
% g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
% g.set_text_options('interpreter','tex','base_size',9);
% g.no_legend();
% g.draw();
% g.redraw(0.025);
% 
clc; clear g; figure('position',[0 0 500 300]);
g = gramm('x',(toy_promoter_response.(var_to_plot)),'y',(toy_promoter_response.mCitrine_max),...
    'color',toy_promoter_response.Msn2_tex,...
    'subset',ismember(toy_promoter_response.pulse_idx,conditions_to_plot) & ...
    toy_promoter_response.param_scaled==param_scaled & toy_promoter_response.n==1 & ...
    ismember(toy_promoter_response.K,K_list) & ismember(toy_promoter_response.scale_factor,scale_factor_list));
g.facet_grid(toy_promoter_response.scale_factor,toy_promoter_response.K,'scale','independent');
% g.geom_line();
g.stat_summary();
g.set_names('x',x_label,'y','','column','','row','','color','')
g.set_order_options('row',row_order,'color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.set_text_options('interpreter','tex','base_size',9);
g.no_legend();
g.draw();
g.redraw(0.05);

clc; clear g; figure('position',[0 0 500 300]);
g = gramm('x',(toy_promoter_response.(var_to_plot)),'y',(toy_promoter_response.mCitrine_max_vs_WT),...
    'color',toy_promoter_response.Msn2_tex,...
    'subset',ismember(toy_promoter_response.pulse_idx,conditions_to_plot) & ...
    toy_promoter_response.param_scaled==param_scaled & toy_promoter_response.n==1 & ...
    ismember(toy_promoter_response.K,K_list) & ismember(toy_promoter_response.scale_factor,scale_factor_list));
g.facet_grid(toy_promoter_response.scale_factor,toy_promoter_response.K,'scale','fixed');
% g.geom_line();
g.stat_summary();
g.set_names('x',x_label,'y','rel. mCitrine','column','','row','','color','')
g.set_order_options('row',row_order,'color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.set_text_options('interpreter','tex','base_size',9);
g.no_legend();
g.draw();
g.redraw(0.025);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
