function [ MeanAbsErr,STDAbsErr ] = ErrWRTSample( data, MaxTACP, delays, plotit ,Errsize)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
%% Section 3.2: Examine magnitude of error WRT Optimal for both conditions    
   
% This section calculates and plots the magnitude of error with respect to
% the target for Per trials and Mean for Pre Trials as a fxn of TACP and Delay.

% MeanAbsErr should remain relativley stable with TACP but increase for delay for Per trials
%   it should increase with delay and decrease with TACP for Pre Trials
% Same should be true for STD AbsErr


%Plot with X axis TACP, useful to see effect of TACP

MeanAbsErr=[];
STDAbsErr=[];

    for type=1:2
        for delay=1:length(delays)
    for TACP=1:MaxTACP        
Errors=abs(data{data{:,2}==type & data.TACP==TACP & data.DelayTime==delays(delay),18});
trueerr=data{data{:,2}==type & data.TACP==TACP & data.DelayTime==delays(delay),18};
Errors=Errors(Errors<Errsize);
trueerr=trueerr(trueerr<Errsize);
MeanAbsErr(TACP,delay,type)=nanmean(Errors);
STDAbsErr(TACP,delay,type)=nanstd(trueerr);

    end
    
    if plotit
    figure(5)
    hold on
        errorbar([1:MaxTACP],MeanAbsErr(:,delay+3*(type-1)),STDAbsErr(:,delay+3*(type-1)));
        legend('Per 1', 'Per 2', 'Per 6', 'Pre 1', 'Pre 2', 'Pre 6');
xlabel('TACP');
ylabel('Absolute Error from Optimal');
    figure(6)
    hold on
       plot([1:MaxTACP],STDAbsErr(:,delay+3*(type-1)));
       
       legend('Per 1', 'Per 2', 'Per 6', 'Pre 1', 'Pre 2', 'Pre 6');
xlabel('TACP');
ylabel('STD of Absolute Error from Optimal');
    end
    
    end
    end
    

% Plot with X axis as Delay, useful to see effect of Delay    
    for type=1:2
        
    for TACP=1:MaxTACP   
    
        if plotit
    figure(7)
    hold on
        errorbar(delays,MeanAbsErr(TACP,[1:length(delays)]+3*(type-1)),STDAbsErr(TACP,[1:length(delays)]+3*(type-1)));
        legend('TACP1', 'TACP2', 'TACP3', 'TACP4', 'TACP5', 'TAPC6');
xlabel('Delay(s)');
ylabel('Absolute Error from Optimal');
    figure(8)
    hold on
       plot(delays,STDAbsErr(TACP,[1:length(delays)]+3*(type-1)));
       legend('TACP1', 'TACP2', 'TACP3', 'TACP4', 'TACP5', 'TAPC6');
xlabel('Delay (s)');
ylabel('STD of Absolute Error from Optimal');
        end
        
    end
    end



end

