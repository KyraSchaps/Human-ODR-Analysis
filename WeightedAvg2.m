function [ ProportionMean, MeanResidual ] = WeightedAvg2( data, MaxTACP, plotit, LargeCutOff, Cutoff, SessCut, endsession,delays,Errorsize )
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
     
colors={'r','r','b','b'};    
symbol={'o','*','o','*'};
    
ProportionMean=[];
MeanResidual=[];
for type=1:2
    for delay=1:length(delays)
        for TACP=1:2
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
        if TACP==1
        TrialsOfInterestANOVA= SessCondition & CPCondition & data.TACP==TACP& abs(data.ErrWRTOptimal)<Errorsize & data{:,2}==type & data.DelayTime==delays(delay);
        Title1={'Change-Point Trials'};
        else
                    TrialsOfInterestANOVA= SessCondition & CPCondition & data.TACP>=TACP & data.TACP<MaxTACP & abs(data.ErrWRTOptimal)<Errorsize & data{:,2}==type & data.DelayTime==delays(delay);
Title1={'Non-change-point trials'};
        end
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

if type==1
Title2={'Percpetion Trials'};
else
    Title2={'Prediction Trials'};
end

figure(1)
if plotit
    subplot(2,4,4*(TACP-1)+delay+2*(type-1));
    hold on
    scatter((Means-Targs),(y-Targs),colors{delay+2*(type-1)},symbol{delay+2*(type-1)});
    if 4*(TACP-1)+delay+2*(type-1)<5
     title([Title2; 'Delay ',num2str(delays(delay))],'Fontsize',16);
    end
    axis([-40,40,-40,40]);
    xs=[-40,40];
    ys=xs*ProportionMean(TACP,delay,type);

plot([0,0],[-40,40],'k','LineWidth',2);
plot([-40,40],[0,0],'k','LineWidth',2);
plot([-40,40],[-40,40],'k','LineWidth',2);
plot(xs,ys,colors{delay+2*(type-1)},'LineWidth',3);
set(gca,'FontSize',20)
if delay+2*(type-1)+4*(TACP-1)==5
xlabel(['GM-SM (Deg)'],'Fontsize',18);
ylabel({'Non-CP Trials';'Response-SM (Deg)'},'Fontsize',18);
elseif delay+2*(type-1)+4*(TACP-1)==1
   ylabel({'CP Trials'},'Fontsize',18);
end
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
