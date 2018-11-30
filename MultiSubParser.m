cd '/Users/KAS/Documents/MATLAB/ODR Task Analysis/PeopleData/MainTaskData'
SubjectsToParse=dir('*Sub*');
subjects={SubjectsToParse.name};

MaxTACP=6; %How many TACP do you want to look at
plotit=0; %Do you want to plot?
% colors={'y','m','r','g','c','b','--y','--m','--r','--g','--c','--b',':y',':m',':r',':g',':c','b:'};
colors={'m','r','c','b'};
pupil=0;
delays=[2,5];
Cutoffsize=0;
cutoff=70;
SessCut=0;
endsession=0;
Errsize=45; %was 40;

tic
recalc=1;
reparse=0;


if recalc
    if reparse
data={};
    else
            subjects=ones(size(data,1),1);

    end
for i=1:length(subjects)
    if reparse
     data{i,1}=subjects{i};  %Subject name
    data{i,2}=ParseSubjectDataPeople2(subjects{i},pupil,plotit);    %The Parsed Data for that subject
        tocalc=data{i,2};
    else
         tocalc=data{i,2}{1,1};
    end
    
    
    data{i,3}={ScorerMatt(tocalc,delays)}; %The scores for the different sessions of that subject, organized by delay*Type
            %ScorerMatt(data{i,2});
                        %WeightedAvg2 will give scatter plots according to
                        %CP and NCP trials
   [PropMean,MeanResidual]=WeightedAvg(tocalc,MaxTACP,1,Cutoffsize,cutoff, SessCut, endsession,delays,Errsize);
    data{i,4}={PropMean}; % The weight given to mean of Guess=x*mean+(1-x)*targ
    data{i,5}={MeanResidual}; %Measure of goodness of fit of the weighted avg
    [ErrWRTOPT,STDWRTOPT]=ErrWRTOptimal(tocalc,MaxTACP,plotit,Cutoffsize,cutoff,SessCut, endsession,delays,Errsize);
    data{i,6}={ErrWRTOPT};    % For each TACP,Delay,Type, what is the avg err wrt optimal
    data{i,7}={STDWRTOPT};     % For each TACP,Delay,Type, what is the STD err wrt optimal
    Errfits=fMinLinearfit(tocalc,MaxTACP,Cutoffsize,cutoff,SessCut, endsession, Errsize,delays);  %***** THIS MIGHT BE BETER AS EXP FIT FOR TACP
    data{i,8}={Errfits}; %First is for Per, 2nd for Pre, fit err as fxn of TACP, and Delay
   % HistoMaker(data{i,2},MaxTACP);
   data{i,9}={PropMeanLinFit(PropMean,MaxTACP,delays)};
    
    if pupil
   [Baseline,MaxVBaseline]=PupilAnalyisis(tocalc,1,delays,Errsize); %1 should be plotit
       [ErrfitsMean, ErrfitsChange]=fMinLinearfitPupilStuff(tocalc,MaxTACP,Cutoffsize,cutoff,SessCut, endsession, Errsize,delays);

   data{i,10}={Baseline};% 
   data{i,11}={MaxVBaseline};
   data{i,12}={ErrfitsMean};
   data{i,13}={ErrfitsChange};
    else 
    [ErrWRTSamp,STDWRTSamp]=ErrWRTSample(tocalc,MaxTACP,delays, plotit,45);

    data{i,10}={ErrWRTSamp};    % For each TACP,Delay,Type, what is the avg err wrt optimal
    data{i,11}={STDWRTSamp};     % For each TACP,Delay,Type, what is the STD err wrt optimal
    
 
    end
%Plots
end
if reparse
if pupil
    data=cell2table(data,'VariableNames',{'SubjectFile','Parsed_Data','Scores', 'Prop_Mean', 'Mean_Resid','ErrWRTOPT','STDWRTOPT','Lin_Errfits', 'PropMean_fits', 'MattMeanPupil', 'MaxVBaseline','Linfit_Pup_Mean', 'Linfit_Pup_Change'});

else

data=cell2table(data,'VariableNames',{'SubjectFile','Parsed_Data','Scores', 'Prop_Mean', 'Mean_Resid','ErrWRTOPT','STDWRTOPT','Lin_Errfits','PropMean_fits','ErrWRTSample','STDWRTSample'});
end
end
end
toc
%% Errwrt targ

legendlables={};
j=1;
clear h

for i=1:length(subjects)

for type=1:2
    for delay=1:length(delays)
%         figure(1)
%         plot([1:MaxTACP],PropMean(:,delay,type));
%         hold on
        figure(i)

        thecolor=colors{(2*(type-1)+(delay))};
        plot([1:MaxTACP],data.ErrWRTSample{i,1}(:,delay,type),thecolor);
        hold on
        title(data.SubjectFile(i));
                xlabel('TACP');
        ylabel('Mean Abs Err from targ');

        if ~isnan(data.ErrWRTOPT{i,1}(:,delay,type))
%         figure(10);
%         h(j)=plot([1:MaxTACP],data.ErrWRTOPT{i,1}(:,delay,type),colors{(3*(type-1)+(delay))+6*(i-1)});
        j=j+1;
        hold on
        xlabel('TACP');
        ylabel('Mean Abs Err from targ');
        end
        
    end 
end

  
end
% 
%legend([h(:)],'KAS Per1','KAS Per3', 'KAS Per6', 'KAS Pre1', 'KAS Pre3', 'KAS Pre6');

%% Err wrt optimal scatter delay
ErrWRTShortPer=[];
ErrWRTLongPer=[];
ErrWRTShortPre=[];
ErrWRTLongPre=[];
 for ACP=1:MaxTACP
for i=1:length(subjects)
   
ErrWRTShortPer(i,ACP)=data.ErrWRTOPT{i,1}(ACP,1,1);
ErrWRTLongPer(i,ACP)=data.ErrWRTOPT{i,1}(ACP,2,1);

ErrWRTShortPre(i,ACP)=data.ErrWRTOPT{i,1}(ACP,1,2);
ErrWRTLongPre(i,ACP)=data.ErrWRTOPT{i,1}(ACP,2,2);
  
end
figure(1)
scatter(ErrWRTShortPer(:,ACP),ErrWRTLongPer(:,ACP),'b','o');
hold on
plot(x,x);
% figure(2)
hold on
scatter(ErrWRTShortPre(:,ACP),ErrWRTLongPre(:,ACP),'c','*');
x=[0:1:20];
plot(x,x);
axis([0,25,0,25]);
title(['TACP=',num2str(ACP)]);
ylabel('Err WRT Target, 5seconds');
xlabel('Err WRT Target, 2seconds');
 end

 
 %%
 %% Scatter Prop mean
 colours={'r','m','y','g','c','b'};
ErrWRTShortPer=[];
ErrWRTLongPer=[];
ErrWRTShortPre=[];
ErrWRTLongPre=[];
 for ACP=1:MaxTACP
for i=1:length(subjects)
   
ErrWRTShortPer(i,ACP)=data.Prop_Mean{i,1}(ACP,1,1);
ErrWRTLongPer(i,ACP)=data.Prop_Mean{i,1}(ACP,2,1);

ErrWRTShortPre(i,ACP)=data.Prop_Mean{i,1}(ACP,1,2);
ErrWRTLongPre(i,ACP)=data.Prop_Mean{i,1}(ACP,2,2);
  
end
figure(1)
scatter(ErrWRTShortPer(:,ACP),ErrWRTLongPer(:,ACP),colours{ACP},'o');
hold on
plot(x,x);
% figure(2)
hold on
scatter(ErrWRTShortPre(:,ACP),ErrWRTLongPre(:,ACP),colours{ACP},'filled');
x=[0:1:20];
plot(x,x);
axis([0,1,0,1]);
title(['TACP=',num2str(ACP)]);
ylabel('Prop Mean, 5seconds');
xlabel('Prop Mean, 2seconds');
 end

  %%
 %% Scatter Prop mean vs TACP
 colours={'r','m','y','g','c','b'};
ErrWRTShortPer=[];
ErrWRTLongPer=[];
ErrWRTShortPre=[];
ErrWRTLongPre=[];
 for ACP=1:MaxTACP
for i=1:length(subjects)
   
ErrWRTShortPer(i,ACP)=data.Prop_Mean{i,1}(ACP,1,1);
ErrWRTLongPer(i,ACP)=data.Prop_Mean{i,1}(ACP,2,1);

ErrWRTShortPre(i,ACP)=data.Prop_Mean{i,1}(ACP,1,2);
ErrWRTLongPre(i,ACP)=data.Prop_Mean{i,1}(ACP,2,2);
  
end
figure(1)
scatter(ACP*ones(1,i),ErrWRTLongPer(:,ACP),colours{ACP},'o');
hold on
scatter(ACP*ones(1,i),ErrWRTShortPer(:,ACP),colours{ACP},'filled');

 figure(2)
hold on
scatter(ACP*ones(1,i),ErrWRTShortPre(:,ACP),colours{ACP},'filled');
scatter(ACP*ones(1,i),ErrWRTLongPre(:,ACP),colours{ACP},'o');


title(['TACP=',num2str(ACP)]);
ylabel('Prop Mean, 5seconds');
xlabel('Prop Mean, 2seconds');
 end

 
 
%% Err wrt optimal
if 1
legendlables={};
j=1;
% clear h

for i=1:length(subjects)

for type=1:2
    for delay=1:length(delays)
%         figure(1)
%         plot([1:MaxTACP],PropMean(:,delay,type));
%         hold on
        figure(12+i)
        thecolor=colors{(2*(type-1)+(delay))};
        plot([1:MaxTACP],data.ErrWRTOPT{i,1}(:,delay,type),thecolor);
        hold on
        title(data.SubjectFile(i));
        xlabel('TACP');
        ylabel('Mean Abs Err from Optimal');
        
        if ~isnan(data.ErrWRTOPT{i,1}(:,delay,type))
%         figure(10);
%         h(j)=plot([1:MaxTACP],data.ErrWRTOPT{i,1}(:,delay,type),colors{(3*(type-1)+(delay))+6*(i-1)});
        j=j+1;
        hold on
        xlabel('TACP');
        ylabel('Mean Abs Err from Optimal');
        end
        
    end 
end

  
end

%legend([h(:)],'KAS Per1','KAS Per3', 'KAS Per6', 'KAS Pre1', 'KAS Pre3', 'KAS Pre6');
%%
clear h
%% Proportion Mean
j=1;
for i=1:length(subjects)
    for type=1:2
        for delay=1:length(delays)
    figure(17+i);     
    hold on;
    if ~isempty(data.Prop_Mean{i,1}(:,delay,type))
                high=min(3,delay);
        thecolor=colors{(2*(type-1)+(delay))};
h(j)=plot([1:MaxTACP],data.Prop_Mean{i,1}(:,delay,type),thecolor);
        title(data.SubjectFile(i));
    end
            
%     figure(6)
    hold on;
%     h(j)=plot([1:MaxTACP],data.Prop_Mean{i,1}(:,delay,type),colors{(3*(type-1)+(delay))+6*(i-1)});
        j=j+1;
               xlabel('TACP');
        ylabel('Proportion Mean');
        end
    end
end
%legend([h(:)],'KAS Per1','KAS Per3', 'KAS Per6', 'KAS Pre1', 'KAS Pre3', 'KAS Pre6');

% %% Err WRT Feedback
% for i=4%:length(subjects)
% ErrWRTFB( data{i,2}, MaxTACP, 1,i );
% end

% %% Learning Rate( between -1 and 1)
% for i=1:length(subjects)
% LearingRateGraph( data{i,2}, MaxTACP, 1,i );
% end
end


% %% full model (boxplot)
% 
% LinearFitsCombo=[];
% figure(21)
% for i=1:length(subjects)
% LinearFitsCombo(:,i)=data.Lin_Errfits{i,1}{:,4};
% 
% end
% 
% LinearFitsCombo=LinearFitsCombo';
% boxplot(LinearFitsCombo,'Labels',{data.Lin_Errfits{1,1}{:,1}'});
% LinearFitsCombo=[];
% subjectnames=[];
% item=[];
% hold on
% for i=1:length(subjects)
% hold on
% LinearFitsCombo=[LinearFitsCombo;data.Lin_Errfits{i,1}{:,4}];
% subjectnames=[subjectnames;repmat(data.SubjectFile(i),7,1)];
% item=[item;[1:7]'];
% end
% 
% gscatter(item,LinearFitsCombo,subjectnames);
% 
% %% Prop mean fit (Boxplot)
% 
% PropMeanFitsCombo=[];
% figure(20)
% 
% for i=1:4%length(subjects)
% PropMeanFitsCombo(i,:)=[data.PropMean_fits{i,1}{:,4}];
% 
% end
% 
% PropMeanFitsCombo=PropMeanFitsCombo;
% boxplot(PropMeanFitsCombo,'Labels',{data.PropMean_fits{1,1}{:,1}});
% PropMeanFitsCombo=[];
% subjectnames=[];
% item=[];
% hold on
% 
% 
% for i=1:5%length(subjects)
% hold on
% PropMeanFitsCombo=[PropMeanFitsCombo;[data.PropMean_fits{i,1}{:,4}]'];
% subjectnames=[subjectnames;repmat(data.SubjectFile(i),7,1)];
% item=[item;[1:7]'];
% end
% 
% gscatter(item,PropMeanFitsCombo,subjectnames);
%% Partial model fitting error


for j=2:3
LinearFitsCombo=[]; 
figure(23+j)
for i=1:length(subjects)
LinearFitsCombo(:,i)=cell2mat(data.Lin_Errfits{i,1}{1:4,j});
end

LinearFitsCombo=LinearFitsCombo';
boxplot(LinearFitsCombo,'Labels',{data.Lin_Errfits{1,1}{[1:3,7],1}'});
LinearFitsCombo=[];
subjectnames=[];
item=[];
hold on
title(['type', num2str(j-1), 'CP cutoff ', num2str(cutoff), ' Trial Cutoff ', num2str(endsession)]);

for i=1:length(subjects)
hold on
LinearFitsCombo=[LinearFitsCombo;cell2mat(data.Lin_Errfits{i,1}{1:4,j})];
subjectnames=[subjectnames;repmat(data.SubjectFile(i),4,1)];
item=[item;[1:4]'];
end

gscatter(item,LinearFitsCombo,subjectnames);
end

%% Prop mean fit by conditon
for j=2:3
PropMeanFitsCombo=[];
figure(21+j)

for i=1:length(subjects)
PropMeanFitsCombo(i,:)=([data.PropMean_fits{i,1}{1:4,j}]);

end

PropMeanFitsCombo=PropMeanFitsCombo;
boxplot(PropMeanFitsCombo,'Labels',{data.PropMean_fits{1,1}{[1:3,7],1}});
PropMeanFitsCombo=[];
subjectnames=[];
item=[];
title(['type', num2str(j-1), 'CP cutoff ', num2str(cutoff),' ' , ' Trial Cutoff ', num2str(endsession)]);
hold on


for i=1:length(subjects)
hold on
PropMeanFitsCombo=[PropMeanFitsCombo;[data.PropMean_fits{i,1}{:,j}]'];
subjectnames=[subjectnames;repmat(data.SubjectFile(i),4,1)];
item=[item;[1:4]'];
end

gscatter(item,PropMeanFitsCombo,subjectnames);
end


%% Baseline Pupil fit by conditon
for j=2:3
Baseline=[];
figure(21+j)

for i=1:length(subjects)
Baseline(i,:)=(cell2mat([data.Linfit_Pup_Mean{i,1}{1:4,j}]));

end


boxplot(Baseline,'Labels',{data.Linfit_Pup_Mean{1,1}{[1:3,7],1}});
Baseline=[];
subjectnames=[];
item=[];
% title(['type', num2str(j-1), 'CP cutoff ', num2str(2),' ' , ' Trial Cutoff ', num2str(2)]);
hold on


for i=1:length(subjects)
hold on
Baseline=[Baseline;cell2mat(data.Linfit_Pup_Mean{i,1}{:,j})];
subjectnames=[subjectnames;repmat(data.SubjectFile(i),4,1)];
item=[item;[1:4]'];
end

gscatter(item,Baseline,subjectnames);
end


%% Change Pupil fit by conditon
for j=2:3
Change=[];
figure(23+j)

for i=1:length(subjects)
Change(i,:)=(cell2mat([data.Linfit_Pup_Change{i,1}{1:4,j}]));

end


boxplot(Change,'Labels',{data.Linfit_Pup_Change{1,1}{[1:3,7],1}});
Change=[];
subjectnames=[];
item=[];
% title(['type', num2str(j-1), 'CP cutoff ', num2str(cutoff),' ' , ' Trial Cutoff ', num2str(endsession)]);
hold on


for i=1:length(subjects)
hold on
Change=[Change;cell2mat(data.Linfit_Pup_Change{i,1}{1:4,j})];
subjectnames=[subjectnames;repmat(data.SubjectFile(i),4,1)];
item=[item;[1:4]'];
end

gscatter(item,Change,subjectnames);
end


%%
colors={'m','r','c','b'}
if pupil
    for i=1:length(subjects)
        for Type=1:2
            for delay=1:length(delays)

        figure(i)
plot([1:6],data.MattMeanPupil{i,1}(:,delay,Type),colors{delay+2*(Type-1)});
hold on

title('Mean Pupil over 2 seconds post Stim');
            end
        end
        legend('Per2', 'Per5', 'Pre2', 'Pre5');

    end
end
if pupil
    for i=1:length(subjects)
        for Type=1:2
            for delay=1:length(delays)


        figure(6+i)
plot([1:6],data.MattMeanPupil{i,1}(:,delay,Type),colors{delay+2*(Type-1)});
hold on
legend('Per2', 'Per5', 'Pre2', 'Pre5');
title('Pupil Change over 2 seconds post Stim');
            end
        end
    end
end