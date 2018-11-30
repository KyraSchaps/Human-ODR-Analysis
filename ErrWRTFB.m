function [ MeanAbsErr,STDAbsErr ] = ErrWRTFB( data, MaxTACP, plotit, subject )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
%% Section 3.2: Examine magnitude of error WRT Optimal for both conditions    
   
% This section calculates and plots the magnitude of error with respect to
% the target for Per trials and Mean for Pre Trials as a fxn of TACP and Delay.

% MeanAbsErr should remain relativley stable with TACP but increase for delay for Per trials
%   it should increase with delay and decrease with TACP for Pre Trials
% Same should be true for STD AbsErr
data=data{1,1};
colors={'y','m','r','g','c','b','--y','--m','--r','--g','--c','--b',':y',':m',':r',':g',':c','b:'};


%Plot with X axis TACP, useful to see effect of TACP
delays=[1,3,5,6];
MeanAbsErr=[];
STDAbsErr=[];

    for type=1:2
        for delay=1:4
    for TACP=1:MaxTACP       
        TrialsOfInt=data{:,2}==type & data.TACP==TACP & data.DelayTime==delays(delay) & data.PercOfTwoSTDCorr>.04;
RawErrs=(degAngDiff(data{TrialsOfInt,7},data{TrialsOfInt,2+type}));
Errors=abs(degAngDiff(data{TrialsOfInt,7},data{TrialsOfInt,2+type}));
Errors=Errors(Errors<40);
MeanAbsErr(TACP,delay,type)=nanmean(Errors);
STDAbsErr(TACP,delay,type)=nanstd(Errors);
rawSTD(TACP,delay,type)=nanstd(RawErrs);

    end
    
    if plotit
    figure(5)
    hold on
      high=min(3,delay);
        errorbar([1:MaxTACP],MeanAbsErr(:,delay+3*(type-1)),STDAbsErr(:,delay+3*(type-1)),colors{(3*(type-1)+(high))+6*(subject-2)});
        legend('Per 1', 'Per 2', 'Per 6', 'Pre 1', 'Pre 2', 'Pre 6');
xlabel('TACP');
ylabel('Absolute Error from FB');
  figure(7)
    hold on
       plot([1:MaxTACP],STDAbsErr(:,delay+3*(type-1)),colors{(3*(type-1)+(high))+6*(subject-2)});
       
       legend('Per 1', 'Per 2', 'Per 6', 'Pre 1', 'Pre 2', 'Pre 6');
xlabel('TACP');
ylabel('STD of Abs Error from FB');
    figure(8)
    hold on
       plot([1:MaxTACP],rawSTD(:,delay+3*(type-1)),colors{(3*(type-1)+(high))+6*(subject-2)});
       
       legend('Per 1', 'Per 2', 'Per 6', 'Pre 1', 'Pre 2', 'Pre 6');
xlabel('TACP');
ylabel('STD of Raw Error from FB');
    end
    
    end
    end
    




end

