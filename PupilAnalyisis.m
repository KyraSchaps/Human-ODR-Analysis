function [ MattMeanPupil, MattChangePupil] = PupilAnalyisis( data,plotit,delays,Errsize )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
%Plot of subject averagered-pupil response over 2 second period as a fxn of SACP (by Delay and TaskType)


colors=['r','m','y','g', 'c', 'b'];

MattMeanPupil=[];
MattChangePupil=[];
    for Type=1:2
        for delay=1:length(delays)
for TACP=1:6
    
    TrialsOfInterest=  data.TACP==TACP & abs(data.ErrWRTOptimal)<Errsize & data.DelayTime==delays(delay) & data{:,2}==Type;
    if plotit && nnz(TrialsOfInterest)>0
    figure(8)
    subplot(1,4,delay+2*(Type-1));
    hold on;
    avgPup=nanmean(cell2mat(data{TrialsOfInterest,26}));
    plot( avgPup(500:3000),colors(TACP));
    title(['Delay', num2str(delays(delay)), ' TrialType ', num2str(Type)]);   
   legend('1 SACP','2','3','4','5','6');
   xlabel('Time(ms): 500=delay start');
    end
    MattMeanPupil(TACP,delay,Type)=nanmedian(cell2mat(data{TrialsOfInterest,29})); %27
    MattChangePupil(TACP,delay,Type)=nanmedian(cell2mat(data{TrialsOfInterest,30})); %28
      figure(3)
  hold on; 
    subplot(1,4,delay+2*(Type-1));
    hold on;
    if sum(TrialsOfInterest)>0 & plotit;
    holding=nanmean(cell2mat(data{TrialsOfInterest,26}));
    Baseline=nanmean(holding(500:510));
      plot(avgPup(500:3000)-Baseline,colors(TACP));
     legend('1 SACP','2','3','4','5','6');  
     title(['baseline Delay', num2str(delays(delay)), ' TrialType ', num2str(Type)]);   
    end

end
if plotit
figure(12)
plot([1:6],MattMeanPupil(:,delay,Type));
hold on
legend('Per2', 'Per5', 'Pre2', 'Pre5');
title('Baseline Pupil');
xlabel('TACPS');

figure(13)
plot([1:6],MattChangePupil(:,delay,Type));
hold on
legend('Per2', 'Per5', 'Pre2', 'Pre5');
title('Max-Baseline=Evoked');
xlabel('TACPS');
end
    end
        end
 


end

