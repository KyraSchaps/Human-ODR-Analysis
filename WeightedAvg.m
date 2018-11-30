function [ ProportionMean, MeanResidual ] = WeightedAvg( data, MaxTACP, plotit, LargeCutOff, Cutoff, SessCut, endsession,delays,Errorsize )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%% Section 2.1: Fit Guess as a weighted avg of Mean and PTarg
     %ProporitonMeans is measure of what percent of the guess is mean (rest
     %is previous target).  Expect this measure to stay relativley stable
     %for perception trials, possibly increasing a little, and for it to
     %increase with TACP a lot for Prediction Trials.
     
     %Mean Residual is a measure of how clustered around the prediction
     %based on Mean and Ptarg the guess is.  Expect this number to increase
     %with delay.  
     % selection arrays for small, large cps
CUTOFF = LargeCutOff;% %45; %was 75
Lcp = [[data.CPMagnitude]>CUTOFF];
     
    
    
ProportionMean=[];
MeanResidual=[];
for type=1:2
    for delay=1:length(delays)
        for TACP=1:MaxTACP
           if Cutoff 
               CPCondition=Lcp;
           else
               CPCondition=ones(size(Lcp));
           end
           
           if endsession
               SessCondition=data.TrialInSess>(max(data.TrialInSess)-SessCut);
           else
               SessCondition=ones(size(Lcp));
        end
        TrialsOfInterestANOVA= SessCondition & CPCondition & data.TACP==TACP& abs(data.ErrWRTOptimal)<Errorsize & data{:,2}==type & data.DelayTime==delays(delay);

y=data.Guess(TrialsOfInterestANOVA);
Means=data.GenMean(TrialsOfInterestANOVA);
Targs=data.TargAngle(TrialsOfInterestANOVA);
if ~isempty(Means)
[ProportionMean(TACP,delay,type),~,test]=lsqlin((Means-Targs),(y-Targs),[],[],[],[],0,1);
MeanResidual(TACP,delay,type)=mean(abs(test));
else
ProportionMean(TACP,delay,type)=nan;
MeanResidual(TACP,delay,type)=nan;
end

if plotit
    figure(4+delay+2*(type-1))
    subplot(MaxTACP,1,TACP);
    hold on
    scatter((Means-Targs),(y-Targs));
end
        end
        
%         
%         if plotit
%         figure(1)
%     hold on
%         plot([1:MaxTACP],ProportionMean(:,delay,type));
%         legend('Per 1', 'Per 2', 'Per 6', 'Pre 1', 'Pre 2', 'Pre 6');
% xlabel('TACP');
% ylabel('Proportion Mean');
%     figure(2)
%     hold on
%        plot([1:MaxTACP],MeanResidual(:,delay,type));
%        
%        legend('Per 1', 'Per 3', 'Per 6', 'Pre 1', 'Pre 3', 'Pre 6');
% xlabel('TACP');
% ylabel('Residulal of Weighted Avg');



        end
%     end
% end
end
end
