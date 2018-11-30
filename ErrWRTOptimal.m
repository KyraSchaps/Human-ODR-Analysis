function [ MeanAbsErr,STDAbsErr ] = ErrWRTOptimal( data, MaxTACP, plotit, LargeCutOff, Cutoff, SessCut, endsession, delays,Errsize)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
%% Section 3.2: Examine magnitude of error WRT Optimal for both conditions    
   
% This section calculates and plots the magnitude of error with respect to
% the target for Per trials and Mean for Pre Trials as a fxn of TACP and Delay.

% MeanAbsErr should remain relativley stable with TACP but increase for delay for Per trials
%   it should increase with delay and decrease with TACP for Pre Trials
% Same should be true for STD AbsErr

CUTOFF = LargeCutOff;% %45; %was 75
Lcp = [[data.CPMagnitude]>CUTOFF];
%Plot with X axis TACP, useful to see effect of TACP
MeanAbsErr=[];
STDAbsErr=[];

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
             
        
Errors=abs( data{ SessCondition & CPCondition & data{:,2}==type & data.TACP==TACP & data.DelayTime==delays(delay),20});
trueerr=data{ SessCondition & CPCondition & data{:,2}==type & data.TACP==TACP & data.DelayTime==delays(delay),20};
Errors=Errors(Errors<Errsize);
trueerr=trueerr(trueerr<Errsize);
MeanAbsErr(TACP,delay,type)=nanmean(Errors);  % trying median not mean 10/8/18
STDAbsErr(TACP,delay,type)=nanstd(trueerr);

    end
        end

    end
    

  



end

