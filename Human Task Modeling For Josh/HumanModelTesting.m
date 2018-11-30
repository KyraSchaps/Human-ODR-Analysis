%This script will :
% Describe the modeles tested here

% 0) Set up parameters for model fitting

% 1) Find the STD of error in the Estimation reports from actual subject
% data

% 2) Generate subect-number matched (but not actual samples seen by
% subjects) model predictions for each of the different possible models of
% delay-modified inference 

% 3) For a single simulated subject with more trials, examine how well fitting
% can recapture parameters of simulation

    %Note, all model fits will have an input of the True Sample as the
    %model cannot know what the subject's Estimate of the sample is when
    %fitting.
    
    %Currently this is set up to fit the H rate given the true other
    %parameters of the model, but other parameters can be fit as well


%% Model Descriptions:

%Estimation:
    %Estimation (i.e. subject response)~N(True Sample,DelaySTD);
    %DelaySTD is extracted from subject data

%Inference/Predictions:
    %Model 0: No effect of Delay
        %Pure Reduced Bayesian Model (RBM) from Nassar paper From Rational regulation of learning dynamics by pupil-linked arousal systems
        %Feed in: True Samples experienced by (simulated) subject
        %Parameters (which can be fixed or free for fitting): H and STD of
             %Gen Distribution
        
    %Model 1: Delay affects quality of evidence used by RBM
        %Pure Reduced Bayesian Model (RBM) 
        %Feed in: Estimation of True Samples
        %Parameters (which can be fixed or free for fitting): H and STD of
             %Gen Distribution
             
    %Model 2: Delay affects fidelity of Inference at each step
        %Modified Reduced Bayesian Model (RBM)=CPPandRelia4
            %At each step Belief(t+1)=N(Belief(t),DelaySTD)
        %Feed in: True Samples experienced by (simulated) subject
        %Parameters (which can be fixed or free for fitting): H and STD of
             %Gen Distribution, DelaySTD (currently matched to estimation)
                         
    %Model 3: Delay affects fidelity of Inference at each step and quality
                %of evidence
        %Modified Reduced Bayesian Model (RBM)=CPPandRelia4
        %Feed in: Estimation of True Samples
        %Parameters (which can be fixed or free for fitting): H and STD of
             %Gen Distribution, DelaySTD (currently matched to estimation)             
        

%% 0)  Set up paremters

%load data:
load('PepleParsed.mat');

% data columns:
% 1). Subject name;
% 2). Raw parsed data
    % Contains individual session and trial data for a given subject.
    % Colummns include variables like Sample, Generative Mean (GM),
    % Generative and Sample STD, Subject response (guess), trial type
    % (Est or Prediction), Delay, TACP, pupil during trial. Also included
    % are some summary stats, e.g. mean pupil change and error of response
    % with respect to various reference points like target or GM.
% 3). Scores=performance/payment
% Any structure with a 8x2x2 size is organized (TACP, Delay, Condition) with
% Delay=[2,5], and Condition=[Estimation, Prediction].  tried to name
% reasonanably
% Any structure with a 7x4 size is the coeeficients of a linear regression 
    
    
numsubs=size(data,1);  %Number of subjects actually used, incase you want to model exactly what would see using the models compared to actual
Hreal=.15;             %The H rate used to generate any fake data, currently set to match experiment
STDGen=15;             %The STD of Gen Distribution used to generate any fake data, currently matches experiment
Hsub=.15;              %The H rate used to generate any fake responses, currently set to match the H used to generate the data;
pointnum=200;          %The trials per condition, currently set to match actual experiment
ErrorCutoff=45;        %Exclusion of trials with errror>ErrorCutoff as in actual subjects these are generally lapses and/or blink errors
clear Beliefs;

%% 1) Find the STD of the error in Estimation reports from actual subjects data to be used in modeling


%Extract the ErrWRTSample for estimation task for both delays (short=2,
%long=5) from each subject for each TACP

%Note, any sample with an error of >ErrorCutoff is excluded from further analysis as
%most of these are lapses and/or blink errors

for ACP=1:6
for i=1:numsubs

STDWRTEST(i,ACP,1)=data.STDWRTSample{i,1}(ACP,1,1);
STDWRTEST(i,ACP,2)=data.STDWRTSample{i,1}(ACP,2,1);
  
end
end

%STD used for modeling is the average STD over subjects of the average over
%TACP for each delay condition
DelaySTDSubject=[mean(mean(STDWRTEST(:,:,1))), mean(mean(STDWRTEST(:,:,2)))]
DelaySTD=repmat(DelaySTDSubject,pointnum,1);   

%% 2) Simulate possible subject responses (number matched to actual number of participants) 

for sub=1:numsubs
%% 2.1 Generate a bunch of fake data;

[GM,Samples,CP]=GenFakeCPData2(pointnum,Hreal,STDGen);

%% 2.2 Estimation:

%Assumpiton: TACP does NOT matter for Estimation, only the last sample and
%the delay

    %Estimation (i.e. subject repsonse)~N(True Sample,DelaySTD);


%Format the actual subject STD from above to generate data
                     
doubleS=[Samples,Samples];                                          
Estimations=randn(size(doubleS)).*DelaySTD+doubleS;

%Find the difference between simulated responses (Estimations) and true
%data points

Err=degAngDiff(Estimations,doubleS);

% Get Mean and STD of Error in Estimation for each TACP (accounting for
% instances where there are no samples for a given TACP)

% MeanErrEST(Tacp,sub)(Short,Long)
for TACP=1:6

if length(nanmean(abs(Err(CP==TACP,:))))==2
MeanErrEST(TACP,sub)={[nanmean(abs(Err(CP==TACP & [abs(Err(:,1))<ErrorCutoff]',1))),nanmean(abs(Err(CP==TACP & [abs(Err(:,2))<ErrorCutoff]',2)))]};
elseif length(nanmean(abs(Err(CP==TACP,:))))==1
    MeanErrEST(TACP,sub)={[nanmean(abs(Err(CP==TACP & [abs(Err(:,1))<ErrorCutoff]',:))),nan]};
else
    MeanErrEST(TACP,sub)=[nan,nan];
end
STDSim(TACP,sub,1)=nanstd(Err(CP==TACP & [abs(Err(:,1))<ErrorCutoff]',1));
STDSim(TACP,sub,2)=nanstd(Err(CP==TACP & [abs(Err(:,2))<ErrorCutoff]',2));
end

%% 2.3 Preidictions (Belief)

%Beliefs(Delay,Trial,Model+1). Delay=[2,5], Model=0-3

%Model 0: 

%Generate Responses, should be the same because delay doesn't matter
Beliefs(1,:,1)=CPPandRelia(Samples',Hsub,STDGen); 
Beliefs(2,:,1)=CPPandRelia(Samples',Hsub,STDGen);


%Error of belief=response and true GM, as examined for actual subjects
[ErrPred]=[degAngDiff(Beliefs(1,:,1),GM');degAngDiff(Beliefs(2,:,1),GM')];

%Extract average error of model belief vs true GM
for TACP=1:6
    if length(nanmean(abs(ErrPred(:,CP==TACP)')))==2
        MeanErrPre0(TACP,sub)={[nanmean(abs(ErrPred(1,CP==TACP & [abs(ErrPred(1,:))<ErrorCutoff])')),nanmean(abs(ErrPred(2,CP==TACP & [abs(ErrPred(2,:))<ErrorCutoff])'))]};
    elseif length(nanmean(abs(ErrPred(:,CP==TACP)')))==1
      MeanErrPre0(TACP,sub)={[nanmean(abs(ErrPred(1,CP==TACP & [abs(ErrPred(1,:))<ErrorCutoff])')),nan]};
    else
        MeanErrPre0(TACP,sub)=[nan,nan];
    end
end



%Model 1: Assumption, uses the Estimations generated above to update belief, but belief is the one generated by delta rule model.  

%******Note: STD of GM used to generate these model predictions is
%increased compared to true in order to more closely match actual subject
%resposnes***************

Beliefs(1,:,2)=CPPandRelia(Estimations(:,1)',Hsub,STDGen+15);  
Beliefs(2,:,2)=CPPandRelia(Estimations(:,2)',Hsub,STDGen+15);

%Error of belief=response and true GM, as examined for actual subjects
[ErrPred]=[degAngDiff(Beliefs(1,:,2),GM');degAngDiff(Beliefs(2,:,2),GM')];


%Extract average error of model belief vs true GM
for TACP=1:6
if length(nanmean(abs(ErrPred(:,CP==TACP)')))==2
MeanErrPre(TACP,sub)={[nanmean(abs(ErrPred(1,CP==TACP & [abs(ErrPred(1,:))<ErrorCutoff])')),nanmean(abs(ErrPred(2,CP==TACP & [abs(ErrPred(2,:))<ErrorCutoff])'))]};
elseif length(nanmean(abs(ErrPred(:,CP==TACP)')))==1
    MeanErrPre(TACP,sub)={[nanmean(abs(ErrPred(1,CP==TACP & [abs(ErrPred(1,:))<ErrorCutoff])')),nan]};
else
    MeanErrPre(TACP,sub)=[nan,nan];
end
end


%Model 2: User uses true SM for each update, but Belief is degraded over the
%delay
Beliefs(1,:,3)=CPPandRelia4(Samples',Hsub,STDGen,DelaySTD(1,1));
Beliefs(2,:,3)=CPPandRelia4(Samples',Hsub,STDGen,DelaySTD(1,2));


[ErrPred]=[degAngDiff(Beliefs(1,:,3),GM');degAngDiff(Beliefs(2,:,3),GM')];

for TACP=1:6
if length(nanmean(abs(ErrPred(:,CP==TACP)')))==2
if length(nanmean(abs(ErrPred(:,CP==TACP)')))==2
MeanErrPre2(TACP,sub)={[nanmean(abs(ErrPred(1,CP==TACP & [abs(ErrPred(1,:))<ErrorCutoff])')),nanmean(abs(ErrPred(2,CP==TACP & [abs(ErrPred(2,:))<ErrorCutoff])'))]};
elseif length(nanmean(abs(ErrPred(:,CP==TACP)')))==1
    MeanErrPre2(TACP,sub)={[nanmean(abs(ErrPred(1,CP==TACP & [abs(ErrPred(1,:))<ErrorCutoff])')),nan]};
else
    MeanErrPre2(TACP,sub)=[nan,nan];
end
end


%Model 3:  Use Post-Delay Estimations and Degrading fidelity prior 
Beliefs(1,:,4)=CPPandRelia4(Estimations(:,1)',Hsub,STDGen,DelaySTD(1,1));
Beliefs(2,:,4)=CPPandRelia4(Estimations(:,2)',Hsub,STDGen,DelaySTD(1,2));


[ErrPred]=[degAngDiff(Beliefs(1,:,4),GM');degAngDiff(Beliefs(2,:,4),GM')];
for TACP=1:6
if length(nanmean(abs(ErrPred(:,CP==TACP)')))==2
MeanErrPre3(TACP,sub)={[nanmean(abs(ErrPred(1,CP==TACP & [abs(ErrPred(1,:))<ErrorCutoff])')),nanmean(abs(ErrPred(2,CP==TACP & [abs(ErrPred(2,:))<ErrorCutoff])'))]};
elseif length(nanmean(abs(ErrPred(:,CP==TACP)')))==1
    MeanErrPre3(TACP,sub)={[nanmean(abs(ErrPred(1,CP==TACP & [abs(ErrPred(1,:))<ErrorCutoff])')),nan]};
else
    MeanErrPre3(TACP,sub)=[nan,nan];
end
end
temp=[MeanErrPre3{:,sub}];

 
end
%Check to see the STD coming back on the Estimates matches that going in
STDsimdelay=[mean(mean(STDSim(:,:,1),2)),mean(mean(STDSim(:,:,2),2))]
end
%% 2.4 Cross Subjects summary statistics and figure generating
 EST=nan(6,numsubs,2);
 Pre0=nan(6,numsubs,2);
 Pre1=nan(6,numsubs,2);
 Pre2=nan(6,numsubs,2);
 Pre3=nan(6,numsubs,2);
   

for TACP=1:6
    tempEST=[MeanErrEST{TACP,:}];
    temp2=tempEST(1,1:2:end);
    EST(TACP,:,1)=tempEST(1,1:2:end);
    EST(TACP,:,2)=tempEST(1,2:2:end);
    
    
    tempEST=[MeanErrPre0{TACP,:}];
    Pre0(TACP,:,1)=tempEST(1,1:2:end);
    Pre0(TACP,:,2)=tempEST(1,2:2:end);
    
    tempEST=[MeanErrPre{TACP,:}];
    Pre1(TACP,:,1)=tempEST(1,1:2:end);
    Pre1(TACP,:,2)=tempEST(1,2:2:end);
    
    
        tempEST=[MeanErrPre2{TACP,:}];
    Pre2(TACP,:,1)=tempEST(1,1:2:end);
    Pre2(TACP,:,2)=tempEST(1,2:2:end);
    
    tempEST=[MeanErrPre3{TACP,:}];
    Pre3(TACP,:,1)=tempEST(1,1:2:end);
    Pre3(TACP,:,2)=tempEST(1,2:2:end);
    
 
    
end

%%
MultiSubMeanErrEST=[nanmean(EST(:,:,1),2),nanmean(EST(:,:,2),2)];
MultiSubstdErrEST=[nanstd(EST(:,:,1)',0,1)'/(numsubs-1),nanstd(EST(:,:,2)',0,1)'/(numsubs-1)];

MultiSubMeanErrPre0=[nanmean(Pre0(:,:,1),2),nanmean(Pre0(:,:,2),2)];
MultiSubstdErrPre0=[nanstd(Pre0(:,:,1)',0,1)'/(numsubs-1),nanstd(Pre0(:,:,2)',0,1)'/(numsubs-1)];

MultiSubMeanErrPre1=[nanmean(Pre1(:,:,1),2),nanmean(Pre1(:,:,2),2)];
MultiSubstdErrPre1=[nanstd(Pre1(:,:,1)',0,1)'/(numsubs-1),nanstd(Pre1(:,:,2)',0,1)'/(numsubs-1)];

MultiSubMeanErrPre2=[nanmean(Pre2(:,:,1),2),nanmean(Pre2(:,:,2),2)];
MultiSubstdErrPre2=[nanstd(Pre2(:,:,1)',0,1)'/(numsubs-1),nanstd(Pre2(:,:,2)',0,1)'/(numsubs-1)];

MultiSubMeanErrPre3=[nanmean(Pre3(:,:,1),2),nanmean(Pre3(:,:,2),2)];
MultiSubstdErrPre3=[nanstd(Pre3(:,:,1)',0,1)'/(numsubs-1),nanstd(Pre3(:,:,2)',0,1)'/(numsubs-1)];
%%

figure(1)

subplot(1,3,1);

errorbar([1:TACP],MultiSubMeanErrEST(:,1)',MultiSubstdErrEST(:,1),':or','LineWidth',2);
hold on
errorbar([1:TACP],MultiSubMeanErrEST(:,2)',MultiSubstdErrEST(:,2),'-*r','LineWidth',2);
errorbar([1:TACP],MultiSubMeanErrPre1(:,1)',MultiSubstdErrPre1(:,1),':ob','LineWidth',2);
errorbar([1:TACP],MultiSubMeanErrPre1(:,2)',MultiSubstdErrPre1(:,2),'-*b','LineWidth',2);

legend('2 sec EST','5 sec EST','2 sec Pre','5 sec Pre');
%axis([0,7,8,17]);
title('Model 1');
ylabel({'Simulated Mean Absolute Err from'; 'SM (Estimation) or'; 'GM (Prediction)'});
set(gca,'FontSize',20)
grid on


subplot(1,3,2);

errorbar([1:TACP],MultiSubMeanErrEST(:,1)',MultiSubstdErrEST(:,1),':or','LineWidth',2);
hold on
errorbar([1:TACP],MultiSubMeanErrEST(:,2)',MultiSubstdErrEST(:,2),'-*r','LineWidth',2);
errorbar([1:TACP],MultiSubMeanErrPre2(:,1)',MultiSubstdErrPre2(:,1),':ob','LineWidth',2);
errorbar([1:TACP],MultiSubMeanErrPre2(:,2)',MultiSubstdErrPre2(:,2),'-*b','LineWidth',2);

legend('2 sec EST','5 sec EST','2 sec Pre','5 sec Pre');
title('Model 2');
%axis([0,7,8,17]);
xlabel('TACP');
set(gca,'FontSize',20)
grid on



subplot(1,3,3);

grid on
errorbar([1:TACP],MultiSubMeanErrEST(:,1)',MultiSubstdErrEST(:,1),':or','LineWidth',2);
hold on
errorbar([1:TACP],MultiSubMeanErrEST(:,2)',MultiSubstdErrEST(:,2),'-*r','LineWidth',2);
errorbar([1:TACP],MultiSubMeanErrPre3(:,1)',MultiSubstdErrPre3(:,1),':ob','LineWidth',2);
errorbar([1:TACP],MultiSubMeanErrPre3(:,2)',MultiSubstdErrPre3(:,2),'-*b','LineWidth',2);
title('Model 3');
legend('2 sec EST','5 sec EST','2 sec Pre','5 sec Pre');
%axis([0,7,8,17]);
set(gca,'FontSize',20)
grid on



%% 3). Single Subject simulation of different models and trying to fit back.
%%Generate a bunch of fake data for one subject with more points;
pointnum=2000;
Hreal=.15;
Hsub=.15;
clear Beliefs;

[GM,Samples,CP]=GenFakeCPData2(pointnum,Hreal,STDGen);
DelaySTD=repmat(DelaySTDSubject,pointnum,1);
doubleS=[Samples,Samples];
%

%Estimations
Estimations=randn(size(doubleS)).*DelaySTD+doubleS;

% Beliefs= (Delay, Trials, Model+1);
%Model 0: Use itterative proces with no degredation in sample or Mean
Beliefs(1,:,1)=CPPandRelia(Samples',Hsub,STDGen); 
Beliefs(2,:,1)=CPPandRelia(Samples',Hsub,STDGen);

%Model 1: Assumption, uses the Estimations generated above to update belief, but belief is the one generated by delta rule model.  

Beliefs(1,:,2)=CPPandRelia(Estimations(:,1)',Hsub,STDGen);
Beliefs(2,:,2)=CPPandRelia(Estimations(:,2)',Hsub,STDGen);

%Model 2: User uses true SM for each update, but Belief is degraded over the
%delay
Beliefs(1,:,3)=CPPandRelia4(Samples',Hsub,STDGen,DelaySTD(1,1));
Beliefs(2,:,3)=CPPandRelia4(Samples',Hsub,STDGen,DelaySTD(1,2));

%Model 3:  Use Post-Delay Estimations and Degrading fidelity prior 
Beliefs(1,:,4)=CPPandRelia4(Estimations(:,1)',Hsub,STDGen,DelaySTD(1,1));
Beliefs(2,:,4)=CPPandRelia4(Estimations(:,2)',Hsub,STDGen,DelaySTD(1,2));




%% 3.1 Which do you want to fit?

%Which data do you want to try fitting?
Delay=2;         %1=2sec, 2=5sec
Model=4;         %Index of which model you want=model#+1

%Which model are you trying to fit to?
FittingTo=1;  %1=no noise in inference process, 2=noise in inference process

y=Beliefs(Delay,:,Model);

if FittingTo==1
UseModel=@(x)CPPandRelia(Samples',x(1),STDGen);
else
UseModel=@(x)CPPandRelia4(Samples',x(1),STDGen,DelaySTD(1,1));
end

%One to fit model-generated choices (y) to model fit choices (fun)
%One to fit process generated GM to model fit choices (fun2);
fun = @(x)sum(degAngDiff(y,UseModel(x)).^2);
fun2 = @(x)sum(degAngDiff(GM,UseModel(x)).^2);

%fmn con
fitParameter=patternsearch(fun,[.05],[],[],[],[],[.05],[1])



 clear errmodel
possH=[0:1:30];
possH=[0.05:.01:1];
for i=1:length(possH)
    errmodel(i)=fun(possH(i));
end
figure(2)
plot(possH,errmodel);
title('Error landscape');
ylabel('error');
xlabel('H-rate');



