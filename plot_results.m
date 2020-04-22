%% Bayesian analysis of COVID-19 pandemic
%  Utilises data provided by from various sources (including Johns-Hopkins University
%   DXY, National Health Commisson (China), Protezione Civile (Italy), 
%   and https://www.corona-data.ch/ (Switzerland)).
%  
% post processing analysis and plotting
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: plotting functions desiged for processing all results from a run on the cluster.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex');
% load population table

%% data store
DATA_DIR = '../Data/'

% load transition matrix for confirmed cases by country
opts = detectImportOptions([DATA_DIR,'COVID19data/data-raw/covid19_sorted_14apr.csv']);
% note: this is a bit verbose, but it is the most effective way I can get
% to ensure the order and number of columns is arbitrary
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
T_d = readtable([DATA_DIR,'COVID19data/data-raw/covid19_sorted_14apr.csv'],opts);

% define prior
prior.num_params = 8; % [\alpha_0,\alpha,\beta,\gamma,\alpha_tau]
prior.p1 = [0.0  ,0.0, 0.0,0.0,0.0,0.0,0.0,0];
prior.p2 = [2.0 ,100.0, 1.0, 1.0, 1.0, 1.0,2,2];
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

%% load output 

lab = {'\alpha_0','\alpha','\beta','\gamma','\delta','\eta','n','\kappa','R_e/R_0','g(C_T)/g(C_0)','\delta/\beta','(\delta+\beta)/\gamma'};
ii = 1;
for j=1:252
    try
        load(['./results_iso3166_R',num2str(j),'/results_smc_intensity_poisson_reg',num2str(j),'.mat'],'results');
        results_array(ii) = results;
        T = length(results_array(ii).data.C);
        C0 = zeros(size(results_array(ii).part_vals,1),1);
        for i=1:length(C0)
            I = find(results_array(ii).part_sim(i,1:T) > 100);
            C0(i) = results_array(ii).part_sim(i,I(1)) + results_array(ii).part_sim(i,T+I(1)) ...
                + results_array(ii).part_sim(i,2*T+I(1));
        end
        g_C0 = results_array(ii).part_vals(:,2).*(1./(1 + ...
            (C0).^(results_array(ii).part_vals(:,7))));
            
         CT = results_array(ii).part_sim(:,T) + results_array(ii).part_sim(:,2*T) ...
                + results_array(ii).part_sim(:,3*T);
         
         g_CT = results_array(ii).part_vals(:,2).*(1./(1 + ...
                (CT).^(results_array(ii).part_vals(:,7))));

        param = results_array(ii).part_vals;   
        R0 = ((param(:,1) + param(:,2).*(1./(1 +(C0).^(param(:,7))))))./((param(:,4) + param(:,3).*param(:,6)));    
        Re = ((param(:,1) + param(:,2).*(1./(1 +(CT).^(param(:,7))))))./((param(:,4) + param(:,3).*param(:,6)));    
        
        results_array(ii).part_vals = [results_array(ii).part_vals,(R0-Re)./R0,(g_C0 - g_CT)./g_C0,param(:,5)./param(:,3),(param(:,3)+ param(:,5))./param(:,4)]
        ii = ii + 1;
    end
end
iso_id = []
for i=1:length(results_array)
    iso_id = [iso_id;results_array(i).ISO3166alpha3]
end
for j=1:length(results_array(1).part_vals(1,:))
    X = [];
    G = [];
    for i=1:length(results_array) 
        X = [X;results_array(i).part_vals(:,j)]; 
        G = [G;repmat(results_array(i).ISO3166alpha3,1000,1)];
    end
    figure(252*j);
    boxplot(X,G,'PlotStyle','compact');
    xlabel('Country Code (ISO-3166 alpha-3)');
    ylabel(['$',lab{j}, ' \sim p(',lab{j},' | \mathcal{D})$'])
    ytickformat('%0.2f')
end

%% point estimates
point_ests = zeros(length(results_array),length(results_array(1).part_vals(1,:)));
CI.low = zeros(length(results_array),length(results_array(1).part_vals(1,:)));
CI.high = zeros(length(results_array),length(results_array(1).part_vals(1,:)));
CI.c95 = zeros(length(results_array),length(results_array(1).part_vals(1,:)));

% compute estimates
for i=1:length(results_array)
    point_ests(i,:) = mean(results_array(i).part_vals);    
    CI.low(i,:) = quantile(results_array(i).part_vals,0.25);
    CI.high(i,:) = quantile(results_array(i).part_vals,0.75);
    CI.c95(i,:) = std(results_array(i).part_vals)*1.96/sqrt(size(results_array(i).part_vals,1));
    
end

%% plot predictions 
fid = 1
J = find(ismember(iso_id, {'USA','ESP','ITA','DEU','CHN','KOR','SGP','TWN','IRN','AUS'}))
iso_id(J)
names = {'Australia','China','Germany','Spain','Iran','Italy','South Korea','Singapore','Taiwan','United States'}

%% plot scatter plots of point estimates
J2 = find(ismember(iso_id,{'CHN','KOR','TWN','SGP'})); % regions to label and colour green
J3 = find(ismember(iso_id,{'USA','ESP','DEU','ITA','CHE','IRN','FRA','GBR','TUR'})); % regions to label and colour red
figure;
subplot(1,3,1)
errorbar(point_ests(:,4),point_ests(:,7),CI.c95(:,7),CI.c95(:,7),CI.c95(:,4),CI.c95(:,4),'o','Color',[1 1 1]/2,'Marker','none','CapSize',0);
hold on;
errorbar(point_ests(J3,4),point_ests(J3,7),CI.c95(J3,7),CI.c95(J3,7),CI.c95(J3,4),CI.c95(J3,4),'or','CapSize',0);
text(point_ests(J3,4),point_ests(J3,7),iso_id(J3));
errorbar(point_ests(J2,4),point_ests(J2,7),CI.c95(J2,7),CI.c95(J2,7),CI.c95(J2,4),CI.c95(J2,4),'og','CapSize',0);
text(point_ests(J2,4),point_ests(J2,7),iso_id(J2));
xlabel('$\gamma$');ylabel('$n$');
xlim([0,1]); ylim([0.4,1.8]);

subplot(1,3,2)
errorbar(point_ests(:,4),point_ests(:,8),CI.c95(:,8),CI.c95(:,8),CI.c95(:,4),CI.c95(:,4),'o','Color',[1 1 1]/2,'Marker','none','CapSize',0);
hold on;
errorbar(point_ests(J3,4),point_ests(J3,8),CI.c95(J3,8),CI.c95(J3,8),CI.c95(J3,4),CI.c95(J3,4),'or','CapSize',0);
text(point_ests(J3,4),point_ests(J3,8),iso_id(J3));
errorbar(point_ests(J2,4),point_ests(J2,8),CI.c95(J2,8),CI.c95(J2,8),CI.c95(J2,4),CI.c95(J2,4),'og','CapSize',0);
text(point_ests(J2,4),point_ests(J2,8),iso_id(J2));
xlabel('$\gamma$');ylabel('$\kappa$');
xlim([0,1]); ylim([0,1.8]);

subplot(1,3,3)
errorbar(point_ests(:,7),point_ests(:,8),CI.c95(:,8),CI.c95(:,8),CI.c95(:,7),CI.c95(:,7),'o','Color',[1 1 1]/2,'Marker','none','CapSize',0);
hold on;
errorbar(point_ests(J3,7),point_ests(J3,8),CI.c95(J3,8),CI.c95(J3,8),CI.c95(J3,7),CI.c95(J3,7),'or','CapSize',0);
errorbar(point_ests(J2,7),point_ests(J2,8),CI.c95(J2,8),CI.c95(J2,8),CI.c95(J2,7),CI.c95(J2,7),'og','CapSize',0);
text(point_ests(J2,7),point_ests(J2,8),iso_id(J2));
text(point_ests(J3,7),point_ests(J3,8),iso_id(J3));
xlabel('$n$');ylabel('$\kappa$');
xlim([0.4,1.8]); ylim([0,1.8]);

figure;
subplot(1,2,1)
errorbar(point_ests(:,12),point_ests(:,11),CI.c95(:,11),CI.c95(:,11),CI.c95(:,12),CI.c95(:,12),'o','Color',[1 1 1]/2,'Marker','none','CapSize',0);
hold on;
errorbar(point_ests(J3,12),point_ests(J3,11),CI.c95(J3,11),CI.c95(J3,11),CI.c95(J3,12),CI.c95(J3,12),'or','CapSize',0);
text(point_ests(J3,12),point_ests(J3,11),iso_id(J3));
errorbar(point_ests(J2,12),point_ests(J2,11),CI.c95(J2,11),CI.c95(J2,11),CI.c95(J2,12),CI.c95(J2,12),'og','CapSize',0);
text(point_ests(J2,12),point_ests(J2,11),iso_id(J2));
xlabel('$(\beta+\delta)/\gamma$');ylabel('$\delta/\beta$');

 subplot(1,2,2)
 errorbar(point_ests(:,9),point_ests(:,10),CI.c95(:,10),CI.c95(:,10),CI.c95(:,9),CI.c95(:,9),'o','Color',[1 1 1]/2,'Marker','none','CapSize',0);
 hold on;
 errorbar(point_ests(J3,9),point_ests(J3,10),CI.c95(J3,10),CI.c95(J3,10),CI.c95(J3,9),CI.c95(J3,9),'or','CapSize',0);
 text(point_ests(J3,9),point_ests(J3,10),iso_id(J3));
 errorbar(point_ests(J2,9),point_ests(J2,10),CI.c95(J2,10),CI.c95(J2,10),CI.c95(J2,9),CI.c95(J2,9),'og','CapSize',0);
 text(point_ests(J2,9),point_ests(J2,10),iso_id(J2));
 xlabel('$(R_0-R_e)/R_0$');ylabel('$[g(C_0)-g(C_T)]/g(C_0)$');

% comment this out to continue and print individual prediction results
return

fid = 10
%J = find(ismember(iso_id, {'USA','ITA','CHN','KOR','TWN','IRN'}))
%J = find(ismember(iso_id, {'ITA'}))
%names = {'China','Iran','Italy','South Korea','Taiwan','United States'}

ii =1;
for j=J'
    results = results_array(j)
    j
    
    %% resampling and forward predictions for unobserved future
    pred_days = 7;
    data_pred = results.data;
    optsf = [];
    optsf.handle = figure(fid+1);
    optsf.line_width = 1;
    optsf.alpha = 0.5;
    optsf.error = 'credint';
    % append 10 extra days for simulation
    data_pred.C = [data_pred.C;zeros(pred_days,1)];
    data_pred.D = [data_pred.D;zeros(pred_days,1)];
    data_pred.R = [data_pred.R;zeros(pred_days,1)];
    T = length(data_pred.C);
    N = size(results.part_vals,1);
    predsims = zeros(2*N,3*T);
    epsilon = zeros(2*N,1);
    for i=1:2*N
        D_s = sim_func(data_pred,results.part_vals(mod(i-1,1000)+1,:));
        predsims(i,:) = smry(D_s);
        
        epsilon(i) = sqrt((results.data.C(end) - D_s.C(length(results.data.C))).^2 ...
                            + (results.data.R(end) - D_s.R(length(results.data.R))).^2 ...
                            +(results.data.D(end) - D_s.D(length(results.data.D))).^2);
    end
    [B,I] = sort(epsilon);
    time_seq = T_d.date(T_d.alpha3 == iso_id(j) & T_d.Province_State == 'total')';
    optsf.x_axis = [time_seq, time_seq(end)+1:time_seq(end)+pred_days];
    optsf.color_area = [128 193 219]./255;    % Blue theme
    optsf.color_line = [ 52 148 186]./255;
    optsf = plot_areaerrorbar(predsims(:,1:T),optsf)
    optsf.color_area = [237,177,32]./255;    % yellow theme
    optsf.color_line = [237,177,32]./(2*255);
    optsf = plot_areaerrorbar(predsims(:,T+1:2*T),optsf)
    optsf.color_area = [243 169 114]./255;    % Orange theme
    optsf.color_line = [236 112  22]./255;
    optsf = plot_areaerrorbar(predsims(:,2*T+1:3*T),optsf)
    hold on
    title(sprintf('%s (%s)',names{ii},results.ISO3166alpha3));
    ii = ii + 1;
    plot(time_seq,results.data.C,'ok');
    plot(time_seq,results.data.R,'sk');
    plot(time_seq,results.data.D,'xk');
    ylabel('number of active cases, recoveries and deaths');
    t = find(results.data.C >= 100);
    xlim([time_seq(t(1)), time_seq(end)+pred_days]);
    box on;
    set(gcf,'renderer','Painters')
    print(sprintf('%s_pred.eps',results.ISO3166alpha3),'-depsc2','-tiff','-r300','-painters')
     fid = fid + 1
    %input('next')
end

