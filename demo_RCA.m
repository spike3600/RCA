% START_mixedModel.m
% this script simulates patterns in two regions in an experiment with
% stimuli that have two different features.
% region1 codes for feature1 (invariant to faeture2) and region2 codes for
% feature2 (invariant to feature1).
% the RDMs are tested against a model that hypothesizes pattern
% dissimilarities to be proportional to differences along both features 1
% and 2.
clear;clc;close all
rng('default')
addpath(genpath('/Users/hnili/Desktop/codes/matlab/rsatoolbox-develop'))

%% control vars
nSubjects=20; % number of subjects simulated
nFeature1 = 4;% number of levels for feature1, e.g. orientation
nFeature2 = 4;% number of levels for feature2, e.g. brightness
nCond = nFeature1*nFeature2; % number of conditions in RDMs
nVoxels = [120 150]; % number of voxels/vertices in ROIs 1 and 2 respectively
neural_noise=0.5; % variance (power) of noise in neural patterns 
correlation_type='pearson'; % type of correlation used

% Specification of clusters of representations for generating patterns of
% activity
clusterSpec{1} = 1;
for i=2:nFeature2+1
    clusterSpec{i} = {neural_noise,nFeature1};
end

%% Defining the model RDMs/RDVs

model_feature1 = 1-kron(eye(4),ones(4)); % model of ROI 1
Condition_matrix=reshape(1:nFeature1*nFeature2,[nFeature1 nFeature2])';
idx = Condition_matrix(:);
model_feature2 = model_feature1(idx,idx); % model of ROI 2
model_comb_rdv = squareform(model_feature1.*model_feature2)'; % intermediate model: multiplication of the two ROIs 
model1_rdv = squareform(model_feature1)';   
model2_rdv = squareform(model_feature2)';
%% Loop over subjects
for subject=[1:nSubjects]
    
    %% Experiment 1: evaluating the effect of measuring connectivity through
    % the lens of a intermediate model
    betas_roi1 = rsa.sim.generateBetaPatterns(clusterSpec, nVoxels(1));
    temp = rsa.sim.generateBetaPatterns(clusterSpec, nVoxels(2));
    % re-order to get patterns with similarity structure concordant with
    % model2
    betas_roi2 = temp(idx,:);

    %% Experiment 1: evaluating model-free and model-based connectivity
    rdv_region1 = pdist(betas_roi1)';
    rdv_region2 = pdist(betas_roi2)';
    Exp1_models_corr(subject)=corr(model1_rdv,model2_rdv,'type',correlation_type); % correlation between the two models
    Exp1_model_free_RC(subject)=corr(rdv_region1,rdv_region2,'type',correlation_type); % correlation between ROI1 and ROI2
    Exp1_ROI1_model1_corr(subject)=corr(rdv_region1,model1_rdv,'type',correlation_type); % correlation between ROI1 and model1
    Exp1_ROI2_model2_corr(subject)=corr(rdv_region2,model2_rdv,'type',correlation_type); % correlation between ROI2 and model2
    Exp1_ROI1_comb_model_corr(subject)=corr(rdv_region1,model_comb_rdv,'type',correlation_type); % correlation between ROI1 and combined model
    Exp1_ROI2_comb_model_corr(subject)=corr(rdv_region2,model_comb_rdv,'type',correlation_type); % correlation between ROI2 and combined model

    %% Experiment 2: evaluating the effect of adding noise (common pattern) to both RoIs as a common source
    % Generating RDMs and adding a common noise component to both ROIs
    nTime=200;   % number of time samples incorporated in this analysis
    common_noise_power=7; % variance (strength) of the noise added to both ROIs
    
    % generating activity patterns and RDMs across the analysis time
    for t=1:nTime
        pats1 = rsa.sim.generateBetaPatterns(clusterSpec, nVoxels(1)); % generating patterns of ROI 1
        temp = rsa.sim.generateBetaPatterns(clusterSpec, nVoxels(2));
        pats2 = temp(idx,:);         % generating patterns of ROI 2

        % Generating the noise pattern to add to the activity patterns of both regions
        % here we generate a noise pattern for ROI 1 and use a transformation matrix
        % (again filled with noise) to geenerate correlated patterns for ROI 2                       
        common_noise1=randn(size(pats1,1),nVoxels(1))*common_noise_power; % Generating a random pattern to add to region 1
        common_noise2=common_noise1*randn(nVoxels); % multiply the random pattern from ROI 1 by a random transformation matrix for ROI 2
        
        % record the noise RDMs to use for plotting
        rdv_noise(:,t)=pdist(common_noise1(:,end));
        
        % add the common noise to the activation patterns of regions
        rdvs_roi1_before(:,t) = pdist(pats1);
        rdvs_roi2_before(:,t) = pdist(pats2); 
        rdvs_roi1_after(:,t) = pdist(pats1+common_noise1);
        rdvs_roi2_after(:,t) = pdist(pats2+common_noise2);        
    end
    
    %% Experiment 2: evaluating model-free and model-based connectivity
    
    % Model-free method (direct correlation between RDVs) across time: it
    % incorrectly shows connectivity under the effect of added common noise
    for time=1:nTime
        exp2_r_before(time,1)=corr(rdvs_roi1_before(:,time),rdvs_roi2_before(:,time),'type',correlation_type);
        exp2_r_after(time,1)=corr(rdvs_roi1_after(:,time),rdvs_roi2_after(:,time),'type',correlation_type);
    end
    Exp2_model_free_RC_before(subject)=mean(exp2_r_before);
    Exp2_model_free_RC_after(subject)=mean(exp2_r_after);
       
    % Model-based method with one model (correlation between time courses after calculating
    % neural-RDM correlations): it incorrectly shows connectivity under the
    % effect of added common noise
    for time=1:nTime
        exp2_r_roi1_before(time,1)=corr(rdvs_roi1_before(:,time),model1_rdv,'type',correlation_type);
        exp2_r_roi2_before(time,1)=corr(rdvs_roi2_before(:,time),model1_rdv,'type',correlation_type);
        exp2_r_roi1_after(time,1)=corr(rdvs_roi1_after(:,time),model1_rdv,'type',correlation_type);
        exp2_r_roi2_after(time,1)=corr(rdvs_roi2_after(:,time),model1_rdv,'type',correlation_type);
    end
    Exp2_1model_based_RC_before(subject)=corr(exp2_r_roi1_before,exp2_r_roi2_before,'type',correlation_type);
    Exp2_1model_based_RC_after(subject)=corr(exp2_r_roi1_after,exp2_r_roi2_after,'type',correlation_type);
    
    % Model-based method with two models (correlation between time courses after calculating
    % neural-RDM correlations): it correctly detects no connectivity as it
    % measures the connectivity based on the hypothesized models rather than direct
    % correlation of RDMs
    for time=1:nTime
        exp2_r_roi1_before(time,1)=corr(rdvs_roi1_before(:,time),model1_rdv,'type',correlation_type);
        exp2_r_roi2_before(time,1)=corr(rdvs_roi2_before(:,time),model2_rdv,'type',correlation_type);
        exp2_r_roi1_after(time,1)=corr(rdvs_roi1_after(:,time),model1_rdv,'type',correlation_type);
        exp2_r_roi2_after(time,1)=corr(rdvs_roi2_after(:,time),model2_rdv,'type',correlation_type);
    end
    Exp2_2models_based_RC_before(subject)=corr(exp2_r_roi1_before,exp2_r_roi2_before,'type',correlation_type);
    Exp2_2models_based_RC_after(subject)=corr(exp2_r_roi1_after,exp2_r_roi2_after,'type',correlation_type);
    
    %Model-based method with two random models for further analyses
        model1_rdv_rand = pdist(rand(16))'; % randomizing the mdoels
        model2_rdv_rand = pdist(rand(16))';
   
    for time=1:nTime
        exp2_r_roi1_before_rand(time,1)=corr(rdvs_roi1_before(:,time),model1_rdv_rand,'type',correlation_type);
        exp2_r_roi2_before_rand(time,1)=corr(rdvs_roi2_before(:,time),model2_rdv_rand,'type',correlation_type);
        exp2_r_roi1_after_rand(time,1)=corr(rdvs_roi1_after(:,time),model1_rdv_rand,'type',correlation_type);
        exp2_r_roi2_after_rand(time,1)=corr(rdvs_roi2_after(:,time),model2_rdv_rand,'type',correlation_type);
    end
    Exp2_2models_based_RC_before_rand(subject)=corr(exp2_r_roi1_before_rand,exp2_r_roi2_before_rand,'type',correlation_type);
    Exp2_2models_based_RC_after_rand(subject)=corr(exp2_r_roi1_after_rand,exp2_r_roi2_after_rand,'type',correlation_type);
    Exp2_2models_random_correlation(subject)=corr(model1_rdv_rand,model2_rdv_rand,'type',correlation_type);
    %% Experiment 3: evaluating the ability of connectivity methods to detect connectivity between
    % a pair of brain areas encoding different representations which are
    % either temporally related (congruent: ROI2 shows information coding after ROI1) or unrelated
    % (incongruent: ROI1 shows information coding after ROI2)
    
    clearvars rdvs_roi1 rdvs_roi2 exp3_r exp3_r_roi1 exp3_r_roi2% clearing activity patterns generated for the
    % above experiment and generating new ones with changing temporal dynamics of information encoding
    
    nTime=200; % Number of time samples incorporated in the analysis
    Congruent_delay_imposed=20; % Number of samples of delay from ROI1 to ROI2 (before jitter)
    Incongruent_delay_imposed=-20; % Number of samples of delay from ROI1 to ROI2 (before jitter): it is negative suggesting that ROI2 starts coding before ROI1
    Analysis_delay=20; % Number of samples of delay between ROI1 and ROI2 in connectivity analysis
    Delay_jitter= 10; % Maximum number of samples jittered when imposing inter-area delay to model subject-subject variability
    jitter=randi([-Delay_jitter Delay_jitter],1); % Jitter added to the connectivity delay to model subject-subject variability
    Coding_time_ROI1=[30:60]; % The time spans when each region encodes the information: should be smaller than "nTime"
    Congruent_coding_time_ROI2=[30:50]+Congruent_delay_imposed+jitter;
    Incongruent_coding_time_ROI2=[30:50]+Incongruent_delay_imposed+jitter; 
    
    % generating representational similarity matrix/vectors of the two regions
    % across time. step 1: first generate a random pattern within each area
    % across the analysis time window
    for t=1:nTime
        rdvs_roi1(:,t) = pdist(randn(nFeature1*nFeature2,nVoxels(1)));
        congruent_rdvs_roi2(:,t) = pdist(randn(nFeature1*nFeature2,nVoxels(2)));
        incongruent_rdvs_roi2(:,t) = pdist(randn(nFeature1*nFeature2,nVoxels(2)));
    end
    
    % step 2: now add the meaningful patterns of ROI 1 to the random
    % patterns of step 1 during the specified encoding span
    for t=Coding_time_ROI1
        pats1 = rsa.sim.generateBetaPatterns(clusterSpec, nVoxels(1));
        rdvs_roi1(:,t) = pdist(pats1);
    end
    
    % step 3: now add the meaningful patterns of ROI 2 to the random
    % patterns of step 1 during the specified encoding span
    for t=Congruent_coding_time_ROI2
        temp = rsa.sim.generateBetaPatterns(clusterSpec, nVoxels(2));
        pats2 = temp(idx,:);
        congruent_rdvs_roi2(:,t) = pdist(pats2);
    end
    % again for the incongruent case
    for t=Incongruent_coding_time_ROI2
        temp = rsa.sim.generateBetaPatterns(clusterSpec, nVoxels(2));
        pats2 = temp(idx,:);
        incongruent_rdvs_roi2(:,t) = pdist(pats2);
    end
    
    % Step 4: smoothing the RDVs
    smoothing_window=10;
    rdvs_roi1=smoothdata(rdvs_roi1,2,'movmean',smoothing_window);
    congruent_rdvs_roi2=smoothdata(congruent_rdvs_roi2,2,'movmean',smoothing_window);
    incongruent_rdvs_roi2=smoothdata(incongruent_rdvs_roi2,2,'movmean',smoothing_window);

    %% Experiment 3: evaluating model-free and model-based connectivity

    % Model-free method (direct correlation between RDVs): it misses the
    % connectivity as the two ROIs encode distinct (even opposite information)   
    for time=1:nTime-Analysis_delay
        exp3_r(time,1)=corr(rdvs_roi1(:,time),congruent_rdvs_roi2(:,time+Analysis_delay),'type',correlation_type);
    end
    Exp3_model_free_congruent_RC(subject)=mean(exp3_r);
    
    % Model-free method (direct correlation between RDVs): it misses the
    % connectivity as the two ROIs encode distinct (even opposite
    % information) and the timing is incongruent
    for time=1:nTime-Analysis_delay
        exp3_r(time,1)=corr(rdvs_roi1(:,time),incongruent_rdvs_roi2(:,time+Analysis_delay),'type',correlation_type);
    end
    Exp3_model_free_incongruent_RC(subject)=mean(exp3_r);
    
    
    % 1-Model-based method (correlation between time courses after calculating
    % RDM-model correlations) with congruent representational dynamics: it does not detect  the connectivity because
    % we use only the model from ROI1 which can not capture codes in ROI 2
    % and also the time courses of information is incongruent across areas
    % (ROI2 (destination) before ROI1 (source))
    for time=1:nTime-Analysis_delay
        exp3_r_roi1(time,1)=corr(rdvs_roi1(:,time),model1_rdv,'type',correlation_type);
        exp3_r_roi2(time,1)=corr(incongruent_rdvs_roi2(:,time+Analysis_delay),model1_rdv,'type',correlation_type);
    end
    Exp3_1model_based_incongruent_RC(subject)=corr(exp3_r_roi1,exp3_r_roi2,'type',correlation_type);
   
    
    % 1-Model-based method (correlation between time courses after calculating
    % RDM-model correlations) with congruent representational dynamics: it does not detect  the connectivity because
    % we use only the model from ROI1 which can not capture codes in ROI 2
    for time=1:nTime-Analysis_delay
        exp3_r_roi1(time,1)=corr(rdvs_roi1(:,time),model1_rdv,'type',correlation_type);
        exp3_r_roi2(time,1)=corr(congruent_rdvs_roi2(:,time+Analysis_delay),model1_rdv,'type',correlation_type);
    end
    Exp3_1model_based_congruent_RC(subject)=corr(exp3_r_roi1,exp3_r_roi2,'type',correlation_type);
    
    
    % 2-Models-based method (correlation between time courses after calculating
    % RDM-model correlations) with congruent representational dynamics: it does not detect  the connectivity because
    % the time courses of information is incongruent across areas
    % (ROI2 (destination) before ROI1 (source))
    for time=1:nTime-Analysis_delay
        exp3_r_roi1(time,1)=corr(rdvs_roi1(:,time),model1_rdv,'type',correlation_type);
        exp3_r_roi2(time,1)=corr(incongruent_rdvs_roi2(:,time+Analysis_delay),model2_rdv,'type',correlation_type);
    end
    Exp3_2models_based_incongruent_RC(subject)=corr(exp3_r_roi1,exp3_r_roi2,'type',correlation_type);
    
    
    % 2-Models-based method (correlation between time courses after calculating
    % RDM-model correlations) with congruent representational dynamics: it detects  the connectivity as the two ROIs
    % encode distinct (even opposite information) but at the same time (or with some delay)
    for time=1:nTime-Analysis_delay
        exp3_r_roi1(time,1)=corr(rdvs_roi1(:,time),model1_rdv,'type',correlation_type);
        exp3_r_roi2(time,1)=corr(congruent_rdvs_roi2(:,time+Analysis_delay),model2_rdv,'type',correlation_type);
    end
    Exp3_2models_based_congruent_RC(subject)=corr(exp3_r_roi1,exp3_r_roi2,'type',correlation_type);
    
    fprintf('Subject\t %i%s\n',subject,'  simulated.')
end

%% Plotting and printing the representations and connectivities
%% Experiment 1: Plotting the model and neural RDMs used throughout this paper
title_font_size=14;
graph_pos_and_size=[0.2 0.2 0.6 0.53];% xpos, ypos, xlength, ylength

figure('Color','w','Name','Exp1: model and neural RDMs');
subplot(221);imagesc(model_feature1);axis square off
title('ROI_1 model RDM','fontsize',title_font_size)

subplot(222);imagesc(model_feature2);axis square off
title('ROI_2 model RDM','fontsize',title_font_size)

subplot(223);imagesc(squareform(pdist(betas_roi1)));axis square off
title('ROI_1 neural RDM','fontsize',title_font_size)

subplot(224);imagesc(squareform(pdist(betas_roi2)));axis square off
title('ROI_2 neural RDM','fontsize',title_font_size)

figure('Color','w','Name','Exp1: intermediate model RDM');
gca = axes('Position',graph_pos_and_size); % xpos, ypos, xlength, ylength
imagesc(squareform(model_comb_rdv))
axis square off
title({'Experiment 1:','intermediate model RDM'},'fontsize',title_font_size)

%% Experiment 2: plotting scatter plots with fit
title_font_size=17;
axes_font_size=17;
graph_pos_and_size=[0.22 0.2 0.75 0.66];% xpos, ypos, xlength, ylength
figure('Color','w','Name','Exp2: effect of correlation of models');
gca = axes('Position',graph_pos_and_size); % xpos, ypos, xlength, ylength

x=[Exp2_2models_based_RC_after_rand-Exp2_2models_based_RC_before_rand]';
y=Exp2_2models_random_correlation';
f=polyfit(x,y,1);
scatter(x,y,50,'marker','o','MarkerFaceColor','k');
% set(gca,'fontsize',axes_font_size);

y_est = polyval(f,x);
hold on
plot(x,y_est,'k--','LineWidth',1.5)
hold off
xlabel('Correlation between models [{\itr}]')
ylabel({'Change in connectivity';'[after-before]'});
set(gca,'FontSize',11,'Linewidth',1.5,'fontsize',axes_font_size,...
    'box','off','tickdir','out');

[r,p]=corr(x,y);
text(.2,-0.1,['\itr = ',num2str(r,'%.2f')],'FontSize',16)
text(.2,-0.15,['\itp = ',num2str(p,'%.2f')],'FontSize',16)
title({'Experiment 2:','effect of correlation between models'},'fontsize',title_font_size)


%% Experiment 3: Plotting RDM-Model correlations across time
% calculating the correaltions between the generated RDMs and the desired
% RDM patterns, which should be encoded in each area: just to verify
title_font_size=17;
axes_font_size=17;
graph_pos_and_size=[0.18 0.2 0.75 0.66];% xpos, ypos, xlength, ylength
figure('Color','w','Name','Exp3: temporal dynamics of representations');
gca = axes('Position',graph_pos_and_size); % xpos, ypos, xlength, ylength
for time=1:nTime
    r1(time,1)=corr(rdvs_roi1(:,time),model1_rdv,'type',correlation_type);
    r2(time,1)=corr(congruent_rdvs_roi2(:,time),model2_rdv,'type',correlation_type);
    r3(time,1)=corr(incongruent_rdvs_roi2(:,time),model2_rdv,'type',correlation_type);
end

plot(r1,'color','k','linewidth',2)
hold on
plot(r2,'color',[.6 .6 .6],'linewidth',2)
plot(r3,'--','linewidth',2,'color',[.3 .3 .3])

Additional_gap=0.075;
line([0 nTime-Congruent_delay_imposed],[0 0],'color','k')
xlim([0 nTime-Congruent_delay_imposed])
ylim([min([r1;r2;r3])-Additional_gap max([r1;r2;r3])+Additional_gap])
set(gca,'Linewidth',1.5,'YMinorTick','on','XMinorTick','on','box','off','fontsize',axes_font_size,'box','off','tickdir','out')

title({'Experiment 3:','temporal dynamics of representations'},'fontsize',title_font_size)
legend('ROI_1-model_1','ROI_2-model_2 cng','ROI_2-model_2 inc','Location','NorthEast','edgecolor','none','fontsize',14)
ylabel('Correlation [{\itr}]')
xlabel('Time [sample]')


%% Printing the average connectivity results over subjects
fprintf('Summary results:\n')
fprintf('Exp1:\t\tmodels correlation is %.2f\n',nanmean(Exp1_models_corr))
fprintf('Exp1:\t\tmodel-free RC is %.2f\n',nanmean(Exp1_model_free_RC))
fprintf('Exp1:\t\tstrength of model1 in ROI1 is %.2f\n',nanmean(Exp1_ROI1_model1_corr))
fprintf('Exp1:\t\tstrength of model2 in ROI2 is %.2f\n',nanmean(Exp1_ROI2_model2_corr))
fprintf('Exp1:\t\tmodel-based RC with combined model in ROI1 is %.2f\n',nanmean(Exp1_ROI1_comb_model_corr))
fprintf('Exp1:\t\tmodel-based RC with combined model in ROI2 is %.2f\n',nanmean(Exp1_ROI2_comb_model_corr))
fprintf('Exp1:\t\tp-value of ROI-model correlations is %.2f\n',signrank(Exp1_ROI1_comb_model_corr,Exp1_ROI2_comb_model_corr))
fprintf('Exp2:\t\tmodel-free RC before noise is %.2f\n',nanmean(Exp2_model_free_RC_before))
fprintf('Exp2:\t\tmodel-free RC after noise is %.2f\n',nanmean(Exp2_model_free_RC_after))
fprintf('Exp2:\t\tone-model-based RC before noise is %.2f\n',nanmean(Exp2_1model_based_RC_before))
fprintf('Exp2:\t\tone-model-based RC after noise is %.2f\n',nanmean(Exp2_1model_based_RC_after))
fprintf('Exp2:\t\ttwo-models-based RC before noise is %.2f\n',nanmean(Exp2_2models_based_RC_before))
fprintf('Exp2:\t\ttwo-models-based RC after noise is %.2f\n',nanmean(Exp2_2models_based_RC_after))
fprintf('Exp3:\t\tcongruent lagged model-free RC is %.2f\n',nanmean(Exp3_model_free_congruent_RC))
fprintf('Exp3:\t\tincongruent lagged model-free RC is %.2f\n',nanmean(Exp3_model_free_incongruent_RC))
fprintf('Exp3:\t\tcongruent one-model-based RC is %.2f\n',nanmean(Exp3_1model_based_incongruent_RC))
fprintf('Exp3:\t\tincongruent one-model-based RC is %.2f\n',nanmean(Exp3_1model_based_congruent_RC))
fprintf('Exp3:\t\tcongruent two-models-based RC is %.2f\n',nanmean(Exp3_2models_based_incongruent_RC))
fprintf('Exp3:\t\tincongruent two-models-based RC is %.2f\n',nanmean(Exp3_2models_based_congruent_RC))


%% Statistical testing and plotting final connectivity figures
significance_threshold=0.001; % level of sigificance
title_font_size=17;
axes_font_size=17;
Model_free_color=[.3 .3 .3];
One_model_based_color=[.5 .5 .5];
Two_models_based_color=[.7 .7 .7];
graph_pos_and_size=[0.25 0.37 0.6 0.49];% xpos, ypos, xlength, ylength


% Data concatanation
Experiment_1_data=[Exp1_model_free_RC;Exp1_ROI1_comb_model_corr;Exp1_ROI2_comb_model_corr];
Experiment_2_data=[Exp2_model_free_RC_before;Exp2_model_free_RC_after;Exp2_1model_based_RC_before;...
    Exp2_1model_based_RC_after;Exp2_2models_based_RC_before;Exp2_2models_based_RC_after];
Experiment_3_data=[Exp3_model_free_incongruent_RC;Exp3_model_free_congruent_RC;Exp3_1model_based_incongruent_RC;...
    Exp3_1model_based_congruent_RC;Exp3_2models_based_incongruent_RC;Exp3_2models_based_congruent_RC];


% Limits for bar graphs
Additional_gap=0.075;
max_std=max(nanstd([Experiment_1_data;Experiment_2_data;Experiment_3_data]'));
maxy=max(nanmean([Experiment_1_data;Experiment_2_data;Experiment_3_data],2))+max_std+Additional_gap;
miny=min(nanmean([Experiment_1_data;Experiment_2_data;Experiment_3_data],2))-max_std-Additional_gap;


% Experiment 1
figure('Color','w','Name','Exp1: results');
gca = axes('Position',graph_pos_and_size); % xpos, ypos, xlength, ylength
bar(1,nanmean(Experiment_1_data(1,:),2),'facecolor',Model_free_color,'edgecolor',[1 1 1]);
hold on;
bar([2 3],nanmean(Experiment_1_data([2 3],:),2),'facecolor',One_model_based_color,'edgecolor',[1 1 1]);
errorbar(nanmean(Experiment_1_data,2),nanstd(Experiment_1_data'),'linestyle','none','marker','none','capsize',0,'color','k','linewidth',3);
ylim([miny maxy]);
ylabel('Correlation [{\itr}]')
set(gca,'Linewidth',1.5,'YMinorTick','on','fontsize',axes_font_size,'xtick',[1;2;3],'xticklabel',{'ROI_1-ROI_2','ROI_1-Model','ROI_2-Model'},'box','off',...
    'tickdir','out','xticklabelrotation',45);
title({'Experiment 1:','use of intermediate model'},'fontsize',title_font_size)

% significance testing
for i = 1:size(Experiment_1_data,1)
    if signrank(Experiment_1_data(i,:),[],'tail','right') < significance_threshold
        plot(i,mean(Experiment_1_data(i,:))+std(Experiment_1_data(i,:))+Additional_gap,'r*','MarkerSize', 10)
    end
end


% Experiment 2
figure('Color','w','Name','Exp2: results');
gca = axes('Position',graph_pos_and_size); % xpos, ypos, xlength, ylength
bar(1,nanmean(Experiment_2_data(1,:),2),'facecolor',Model_free_color,'edgecolor',[1 1 1]);
hold on;
bar(2,nanmean(Experiment_2_data(2,:),2),'facecolor',Model_free_color,'edgecolor',[1 1 1]);
bar(3,nanmean(Experiment_2_data(3,:),2),'facecolor',One_model_based_color,'edgecolor',[1 1 1]);
bar(4,nanmean(Experiment_2_data(4,:),2),'facecolor',One_model_based_color,'edgecolor',[1 1 1]);
bar(5,nanmean(Experiment_2_data(5,:),2),'facecolor',Two_models_based_color,'edgecolor',[1 1 1]);
bar(6,nanmean(Experiment_2_data(6,:),2),'facecolor',Two_models_based_color,'edgecolor',[1 1 1]);
errorbar(nanmean(Experiment_2_data,2),nanstd(Experiment_2_data'),'linestyle','none','marker','none','capsize',0,'color','k','linewidth',3);
ylim([miny maxy]);
ylabel('Correlation [{\itr}]')
set(gca,'Linewidth',1.5,'YMinorTick','on','fontsize',axes_font_size,'xtick',[1:6]','xticklabel',...
    {'Model-free before','Model-free after','1-model before','1-model after','2-models before','2-models after'},'box','off',...
    'tickdir','out','xticklabelrotation',45)
title({'Experiment 2:','addition of common input'},'fontsize',title_font_size)
% significance testing
for i = 1:size(Experiment_2_data,1)
    if signrank(Experiment_2_data(i,:),[],'tail','right') < significance_threshold
        plot(i,mean(Experiment_2_data(i,:))+std(Experiment_2_data(i,:))+Additional_gap,'r*','MarkerSize', 10)
    end
end


% Experiment 3
figure('Color','w','Name','Exp3: results');
gca = axes('Position',graph_pos_and_size); % xpos, ypos, xlength, ylength
bar(1,nanmean(Experiment_3_data(1,:),2),'facecolor',Model_free_color,'edgecolor',[1 1 1]);
hold on;
bar(2,nanmean(Experiment_3_data(2,:),2),'facecolor',Model_free_color,'edgecolor',[1 1 1]);
bar(3,nanmean(Experiment_3_data(3,:),2),'facecolor',One_model_based_color,'edgecolor',[1 1 1]);
bar(4,nanmean(Experiment_3_data(4,:),2),'facecolor',One_model_based_color,'edgecolor',[1 1 1]);
bar(5,nanmean(Experiment_3_data(5,:),2),'facecolor',Two_models_based_color,'edgecolor',[1 1 1]);
bar(6,nanmean(Experiment_3_data(6,:),2),'facecolor',Two_models_based_color,'edgecolor',[1 1 1]);
errorbar(nanmean(Experiment_3_data,2),nanstd(Experiment_3_data'),'linestyle','none','marker','none','capsize',0,'color','k','linewidth',3);
ylim([miny maxy]);
ylabel('Correlation [{\itr}]')
set(gca,'Linewidth',1.5,'YMinorTick','on','fontsize',axes_font_size,'xtick',[1:6]',...
    'xticklabel',{'Inc model-free','Cng model-free','Inc 1-model','Cng 1-model','Inc 2-models','Cng 2-models'},...
    'box','off','tickdir','out','xticklabelrotation',45)
title({'Experiment 3:','transformed representations'},'fontsize',title_font_size)
% significance testing
for i = 1:size(Experiment_3_data,1)
    if signrank(Experiment_3_data(i,:),[],'tail','right') < significance_threshold
        plot(i,mean(Experiment_3_data(i,:))+std(Experiment_3_data(i,:))+Additional_gap,'r*','MarkerSize', 10)
    end
end
