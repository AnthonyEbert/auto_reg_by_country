%% Bayesian analysis of COVID-19 pandemic
%  Utilises data provided by from various sources (including Johns-Hopkins University
%   DXY, National Health Commisson (China), Protezione Civile (Italy), 
%   and https://www.corona-data.ch/ (Switzerland)).
%  
% Simulation based model extending a stochastic SIR model with latent infectous
% populatio and regulatory mechnisms.
%
%  Inference of model parameters is perfromed using Approximate Bayesian Computation
%  Within the Sequential Monte Carlo framework of Drovandi and Pettit (2011).
%  
% Authors:
%     Christopher Drovandi (c.drovandi@qut.edu.au)
%           School of Mathematical Sciences
%           Science and Engineering Faculty 
%           Queensland University of Technology
%
%     David J. Warne (david.warne@qut.edu.au)
%           School of Mathematical Sciences
%           Science and Engineering Faculty 
%           Queensland University of Technology
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: To be used on a PBS Cluster (or similat) using the bash script 
% submit_all.sh due to the differences in configurations across HPC clusters,
% it is likely that modifications to this *.m file and the shell script will
% need to be performed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Algorithm Initialisation
%
% Initialise Random number generator for reproducibility
%clear global
% country_id defined within the current workspace to select the region to perform
% inference on. The country id is provide by the submit script.
rng(mod(1337,country_id),'twister')

%%
%region_id = 'Italy'; % 
iso_id = '';
province_id = 'total';

%% data store
DATA_DIR = '../Data/'

% load transition matrix for confirmed cases by country
opts = detectImportOptions([DATA_DIR,'COVID19data/data-raw/covid19_sorted.csv']);
for i=1:length(opts.VariableNames)
    switch opts.VariableNames{i}
        case 'alpha3'
            opts.VariableTypes{i} = 'categorical';
        case 'Country_Region'
            opts.VariableTypes{i} = 'categorical';
        case 'Province_State'
            opts.VariableTypes{i} = 'categorical';
        case 'date'
            opts.VariableTypes{i} = 'datetime';
        otherwise
            opts.VariableTypes{i} = 'double';
    end
end
T_d = readtable([DATA_DIR,'COVID19data/data-raw/covid19_sorted.csv'],opts);

% load population table
opts = detectImportOptions([DATA_DIR,'COVID-19_ISO-3166/full_list.csv']);
for i=1:length(opts.VariableNames)
    switch opts.VariableNames{i}
        case 'alpha3'
            opts.VariableTypes{i} = 'categorical';
        case 'Province_State'
            opts.VariableTypes{i} = 'categorical';
        case 'population'
            opts.VariableTypes{i} = 'double';
    end
end
P_d = readtable([DATA_DIR,'COVID-19_ISO-3166/full_list.csv'],opts);
iso_list = unique(P_d.alpha3);
if strcmp(iso_id,'')
    iso_id = iso_list(country_id);
end

%% extract region of interest
I = T_d.alpha3 == iso_id & T_d.Province_State == province_id;
time_seq = T_d.date(I);
Data.D = T_d.deaths(I);
Data.R = T_d.recovered(I);
Data.C = T_d.active(I);
J = P_d.alpha3 == iso_id & P_d.Province_State == province_id;
Data.P = P_d.population(J);

% check simulation is viable
if Data.P <= 0 || sum(Data.C > 0) == 0 || P_d.alpha3(J) == 'cruise'
    % unpopulated region or no cases of COVID-19 reported
    return
end

%% set-up ABC-SMC 
% the number of particles
N = 1000; 
% target ABC tolerance. (If zero, use acceptance probability stopping criteria)
epsilon_final = 0; 
% tuning parameters for ABC-SMC -- set to good initial defaults
a = 0.75; 
c = 0.01;
% minimum acceptance probability used as stopping criteria if espilon_final = 0 
% if p_acc_min = 0 then use epsilon_final
p_acc_min = 0.005;

% define prior
prior.num_params = 8; % [\alpha_0,\alpha,\beta,\gamma,\alpha_tau]
prior.p1 = [0.0  ,0.0, 0.0,0.0,0.0,0.0,0.0,0];
prior.p2 = [2.0 ,100.0, 1.0, 1.0, 1.0, 1.0,2,2];
prior.sampler = @() [unifrnd(prior.p1,prior.p2)]; 
prior.pdf = @(theta) prod([unifpdf(theta,prior.p1,prior.p2)]);
prior.trans_f = @(theta) [theta];
prior.trans_finv = @(theta) [theta];

% user functions

% stochastic simulation (use @simuldata_reg_fD for death tally only version)
sim_func = @simuldata_reg;
%discrepancy metric
dist_func = @(S_d,S_s) sqrt(sum((S_d(:) - S_s(:)).^2));
% summary statistic
smry_func = @smry;

% run SMC sampler
[part_vals, part_sim, part_s, sims,epsilon_t,p_acc_t] = smc_abc_rw(Data,...
                                       sim_func,dist_func,smry_func,prior,N,...
                                       epsilon_final,a,c,p_acc_min);
%% save results
results.part_vals = part_vals;
results.part_sim = part_sim;
results.part_s = part_s;
results.sims = sims;
results.epsilon = epsilon_t;
results.p_acc_t = p_acc_t;
results.data = Data;
results.id = country_id;
results.name = sprintf('%s',province_id);
results.ISO3166alpha3 = iso_id;
% save output
save(['results_smc_intensity_poisson_reg',num2str(results.id),'.mat'],'results');
 
