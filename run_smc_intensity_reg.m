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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Algorithm Initialisation
%
% Initialise Random number generator for reproducibility
%clear global
% country_id defined within the current workspace to select the region to perform
% inference on.
rng(mod(1337,country_id),'twister')

% data store
DATA_DIR = '../Data/'

% load transition matrix for confirmed cases by country
T_d = import_covid19data([DATA_DIR,'COVID19data/data-raw/covid19_sorted.csv']);
JHU_ISO3166 = import_JHUtoA3([DATA_DIR,'COVID-19_ISO-3166/JohnsHopkins-to-A3.csv']);
P_d = import_country_full_list([DATA_DIR,'COVID-19_ISO-3166/full_list.csv']); 

% extract table for Country/Region (e.g., define country_id = 44 for china, or 114 for Italy)
T_country_id = T_d(ismember(T_d.CountryRegion,...
                            JHU_ISO3166.CountryRegion(...
                            JHU_ISO3166.alpha3 == P_d.alpha3(country_id))) & ...
                            T_d.ProvinceState == 'total',:);



%% extract region of interest
Data.C = T_country_id.deaths;
Data.R = T_country_id.recovered;
Data.D = T_country_id.active;
Data.P = P_d.population(country_id);

% check simulation is viable
if Data.P <= 0 || sum(Data.C > 0) == 0 || P_d.alpha3(country_id) == 'cruise'
    % unpopulated region or no cases of COVID-19 reported
    return
end

%% Synthetic data for identifiability test
%Data = simuldata_reg(Data,[0.01,41,0.1,0.05,0.1,1,1/2,1])


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
prior.p2 = [1.0 ,100, 1.0, 1.0, 1.0, 1.0,1,2];
prior.sampler = @() [unifrnd(prior.p1,prior.p2)]; 
prior.pdf = @(theta) prod([unifpdf(theta,prior.p1,prior.p2)]);
prior.trans_f = @(theta) [theta];
prior.trans_finv = @(theta) [theta];
% user functions

% stochastic simulation
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
results.name = char(P_d.CountryRegion(country_id))
results.ISO3166alpha3 = char(P_d.alpha3(country_id))
% save output
save(['results_smc_intensity_poisson_reg',num2str(country_id),'_',results.ISO3166alpha3,'.mat'],'results');
%% load output 
load(['results_smc_intensity_poisson_reg',num2str(country_id),'_',results.ISO3166alpha3,'.mat'],'results');

% plot marginal posterior densities
figure;
lab = {'\alpha_0','\alpha','\beta','\gamma','\delta','\eta','n','\kappa'};
for i = 1:8
    subplot(2,4,i);
    ksdensity(results.part_vals(:,i));
    xlim([prior.p1(i),prior.p2(i)])
    xlabel(lab{i})
end
title(sprintf('%d %s %s',country_id,results.name,results.ISO3166alpha3))

% plot samples of posterior predicitve distribution agaist data
figure;
T = length(results.data.C);
for i=1:1000
    ha = plot((results.part_sim(i,1:T))','-b'); ha.Color(4) = 0.05;
    hold on;
    hb = plot(results.part_sim(i,T+1:2*T)','-','Color',[237,177,32]/255); hb.Color(4) = 0.05;
    hc = plot(results.part_sim(i,2*T+1:3*T)','-r'); hc.Color(4) = 0.05;
end
title(sprintf('%d %s %s',country_id,results.name,results.ISO3166alpha3))
plot(results.data.C,'+k');
plot(results.data.R,'+k');
plot(results.data.D,'+k');
xlabel('time');
ylabel('counts')

