function [ MeanAbsErr,STDAbsErr ] = LearingRateGraph( data, MaxTACP, plotit, subject )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
%% Section 3.2: Examine magnitude of error WRT Optimal for both conditions    
   
% This section calculates and plots the magnitude of error with respect to
% the target for Per trials and Mean for Pre Trials as a fxn of TACP and Delay.

% MeanAbsErr should remain relativley stable with TACP but increase for delay for Per trials
%   it should increase with delay and decrease with TACP for Pre Trials
% Same should be true for STD AbsErr

data=data{1,1};

%Plot with X axis TACP, useful to see effect of TACP
delays=[1,3,6];
MeanLR=[];
STDLR=[];
colors={'y','m','r','g','c','b','--y','--m','--r','--g','--c','--b',':y',':m',':r',':g',':c','b:'};


    for type=1:2
        for delay=1:3
    for TACP=1:MaxTACP       
        TrialsOfInt=data{:,2}==type & data.TACP==TACP & data.DelayTime==delays(delay);
        MeanLR(TACP,delay,type)=nanmean(data.LR(abs(data.LR)<1 & TrialsOfInt));
        STDLR(TACP,delay,type)=nanstd(data.LR(abs(data.LR)<1 & TrialsOfInt));
% 
%           if plotit
%               figure(6)
%               hold on
%             scatter(abs(data.PredictionError(abs(data.LR)<1 & TrialsOfInt)),data.LR(abs(data.LR)<1 & TrialsOfInt));
%           end

        
    end
    
    if plotit
    figure(9)
    hold on
        errorbar([1:MaxTACP],MeanLR(:,delay+3*(type-1)),STDLR(:,delay+3*(type-1)),colors{(3*(type-1)+(delay))+6*(subject-1)});
        legend('Per 1', 'Per 2', 'Per 6', 'Pre 1', 'Pre 2', 'Pre 6');
xlabel('TACP');
ylabel('MeanLR');
   
    end
    
    end
    end
    




end

