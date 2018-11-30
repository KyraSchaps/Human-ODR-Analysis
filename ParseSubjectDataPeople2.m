function [ data ] = ParseSubjectDataPeople2( subject, pupil,plotit)
%% Section 1: Parse Data from all Subjects in "MainTaskData"

%Jobs: 
%   1). Parse data from individual subjects into single object
%   2). Calculate Session Score 
%   3). Clean the eye data and parse out relevent sections

%Makes a single table called "data" and a single structure called "DATA"
%which contin the relevent behavioral and eye info respectively for all subjects.  Some of
%the relevant eye is copied into data.  "DATA" is made using the function
%"TreatData" which uses "Clean Blinks"

%The score for the session is also calculated here in "SCORES".  The score is the avg
%distance from the objective as a percentage of 2 std of the generative
%mean.

% statusData fields:
%     isSampled               0, 1 for probe
%     sampleType              1=perception, 2=prediction
%     samplesThisRun          
%     sampleRunLength
%     sampleIndex
%     pSample                 
%     trialList               
%     nextSampleTAC
%     ITI                     constant ~10 ms
%     durationTarget          how long showing blue dots per trial
%     durationDelay           delay just before sample ~2 sec
%     distMean                current generative mean
%     targetAngle             current target mean
%     guessAngle              current inference
%     distSTD                 generative distribution STD
%     miniSTD                 target (blue dots) STD
%     sinceChangePT           num trials since last CP


% Output data matrix columns are:
%  1     ...     subject-session name
%  2     ...     trial type: 1=perception, 2=prediction
%  3     ...     current target
%  4     ...     next target (for prediction error)
%  5     ...     current generative mean
%  6     ...     current generative STD
%  7     ...     current guess
%  8     ...     trials since last change-point
%  9     ...     magnitude of change point    
% 10     ...     Timestamp of delay start via eyelink, i.e. when the target
%                   goes off
% 11     ...     average of samples in current change point
% 12     ...     Model Based Belief
% 13     ...     Model Based CPP, inversely proportional to Relevance
% 14     ...     Model Based Reliability. 
% 15     ...     Delay Time
% 16     ...     Trial Score (how close to objective, as percent of 2 std of gen  mean)
% 17     ...     Err WRT Target
% 18     ...     Err WRT Mean
% 19     ...     Err WRT Optimal (Target for Per trials, Mean for Pre trials)
% 20     ...     Err WRT Objective (Target for Per trials, Next Targ for Pre trials)
% 21     ...     CurrentTarget-Mean
% 22     ...     Array of pupil diameter from -500ms targ offset to 2000 ms
%                   after targ offset
% 23     ...     Mean pupil size during 0:2000ms around Targ offset
% 24     ...     Change in pupil from avg in first second after targ offset
%                   to second second after targ offset
                         
% 25     ...     The distance of the avg response to an optimal target at same location for
%                   the given trial from teh actual target location
% 26     ...     The optimal-corrected error with respect to optimal, i.e.
%                   the residual error with respect to optimal after any
%                   systematic error as been removed
% 27     ...     The optimal-corrected Targ-Mean i.e.
%                   the difference between the avg response to a given
%                   target location and the mean.


%close all
%clear all
cd '/Users/KAS/Documents/MATLAB/ODR Task Analysis/PeopleData';

reclean=1;
if reclean==1
    
tic
s=dir(fullfile('MainTaskData',subject, '*.mat'));
files={s.name};
data = [];


if pupil
DATA=TreatData(subject,plotit);
end
Scores=cell(length(files),2);

for jj=1:length(files);  %For all the files
    if isempty(strfind(files{jj},'WM ODR'))
    load(fullfile('MainTaskData', subject, files{jj}));
    
    Scores{jj,2}=mean([statusData([statusData.correct]>.02).correct]');
    Scores{jj,1}=files{jj};
    
    len = length([statusData.isGoodTrial]);
    if len > 20
                holding2=[statusData.sinceChangePT]';

        

        currentmean=[statusData.targGenDistMean]';
        
        subject_data = [...
            jj.*ones(len, 1), ...
            [statusData.sampleType]', ...
            [statusData.targetAngle]', ...
            [statusData.NextAngle]', ...
            currentmean, ...
            [statusData.distSTD]', ...
            [statusData.guessAngle]', ...
            holding2, ...
            NaN(len, 1),...
            [statusData.eyelinkTSDelay]', ...
            NaN(len, 1),...
            NaN(len, 1),...
            NaN(len, 1),...
            NaN(len, 1),...
            [statusData.durationDelay]'...
            [statusData.correct]'...
            [1:size([statusData.correct]')]'
            ]; 
        
        % compute & save magnitude of previous change-point
        cp_trials = [find(subject_data(:,8)==1); len+1];
        subject_data(1:cp_trials(2)-1,9) = 0; % first values
        for ff = 2:length(cp_trials)-1
            subject_data(cp_trials(ff):(cp_trials(ff+1)-1),9) = ...
                abs(rad2deg(angdiff( ...
                deg2rad(subject_data(cp_trials(ff),4)), ...
                deg2rad(subject_data(cp_trials(ff-1),4)))));
        end
%         
      for i=1:length([statusData.isGoodTrial])
           subject_data(i,11)=meanangle(statusData(i).targPerCP); 
      end
%            disp(i);
               ProbeSamples=[statusData.targetAngle];  % If only want samples since last Probe Trial: statusData(i).targPerSample;
           [Beleif,Relevance,Reliability]=CPPandRelia2(ProbeSamples,statusData(1).H,statusData(1).distSTD);
           subject_data(:,12)=Beleif;
           subject_data(:,13)=Relevance;
           subject_data(:,14)=Reliability;

%        
        
        % add to full data matrix
        data = cat(1, data, subject_data);
        
        % increment subject index
 
    end
end
end
Interm=num2cell(data);
%
for i=1:length(files)
Interm((data(:,1)==i),1)={files(i)};
end


data=cell2table(Interm,'VariableNames',{'Session','TrialType_PercOrPred',...
                                             'TargAngle', 'NextAngle','GenMean', 'GenSTD', ...
                                             'Guess', 'TACP', 'CPMagnitude', 'Eyleink_TimeStamp', 'CP_Sample_Avg', ...
                                             'Model_Belief', 'Model_CPP', 'Model_Reliability', 'DelayTime', 'PercOfTwoSTDCorr', 'TrialInSess'...
                                             });




end


% Make some new data catagories about Error vs Various locations
Means=data.GenMean;

data.ErrWRTTarg=degAngDiff(data.TargAngle,data.Guess);
data.ErrWRTMean=degAngDiff(Means,data.Guess);

WRT=zeros(length(data.ErrWRTTarg),1);
WRT(data.TrialType_PercOrPred==1)=1;
WRT(data.TrialType_PercOrPred==2)=2;
GenM=data.GenMean(WRT==2);
data.ErrWRTOptimal=zeros(length(data.ErrWRTTarg),1);
data.ErrWRTOptimal(WRT==2)=degAngDiff(GenM,data.Guess(WRT==2));
data.ErrWRTOptimal(WRT==1)=degAngDiff(data.TargAngle(WRT==1),data.Guess(WRT==1));

data.ErrWRTObjective=zeros(length(data.ErrWRTTarg),1);
data.ErrWRTObjective(WRT==2)=degAngDiff(data.NextAngle(WRT==2),data.Guess(WRT==2));
data.ErrWRTObjective(WRT==1)=degAngDiff(data.TargAngle(WRT==1),data.Guess(WRT==1));

PE=degAngDiff(data.NextAngle,data.Guess);
data.PredictionError=[nan;PE(1:end-1)];
holding=nan;
for i=2:length(data.Guess)
    holding(i-1)=degAngDiff(data.Guess(i),data.Guess(i-1));
end
data.GuessChange=[nan;holding'];
data.LR=data.GuessChange./data.PredictionError;

%Make all angle points relative to target, e.g. X=-6 if X is 354 and
%targ is 3
data.MeanTargDiff=[degAngDiff(data.TargAngle,data.GenMean)];
data.Guess=data.TargAngle+data.ErrWRTTarg;
data.GenMean=data.TargAngle+data.MeanTargDiff;



if pupil
   %   ************Add the Eye Processing Info
       num_subjects = size(unique(data{:,1}),1);
       fileNAMES=unique(data{:,1});

AllPupilDataOfInt={};
for subject2=1:num_subjects
    ToZeeScore=[];
    Testing=[];
    forColumnAvg=[];
           DATA(subject2).DelayPupilData={};
        
  
            DATA(subject2).DelayStartTimes= data{strcmp(data{:,1},fileNAMES{subject2}) ,10}'; 
            
            [C,B,DelayIndexes]=intersect(DATA(subject2).DelayStartTimes,DATA(subject2).time);
          
            for j=1:length(DelayIndexes)
                if nnz(isnan(DATA(subject2).missingSamp(DelayIndexes(j)-1000:DelayIndexes(j)+3000)))<.3*4000
                DATA(subject2).DelayPupilData{j,1}=[DATA(subject2).zCPupil(DelayIndexes(j)-1000:DelayIndexes(j)+3000)];
                else
                  DATA(subject2).DelayPupilData{j,1}=nans(4001,1);

                end
                ToZeeScore=[ToZeeScore,DATA(subject2).DelayPupilData{j,1}];
            end
        % Zscore accross all the fixation times of interest    
               meanis=nanmean(ToZeeScore);
                stdis=nanstd(ToZeeScore);
            for j=1:length(DelayIndexes)
         
                Testing=(ToZeeScore(:,j)-meanis(j))/stdis(j);
               DATA(subject2).DelayPupilData{j,1}=[Testing]';

%                 DATA(subject2).DelayPupilData{j,1}=[Testing(1+(j-1)*4000:4001+(j-1)*4000)];
%                 forColumnAvg=[forColumnAvg;Testing(1+(j-1)*4000:4001+(j-1)*4000)];
            end
%        %Subtract the mean of the column to account for within trial effects     
%         Testing=nanmean(forColumnAvg,1);
%         for i=1:length(Testing)
%             forColumnAvg(:,i)=forColumnAvg(:,i)-Testing(i);
%         end
%         
%             for j=1:length(DelayIndexes)
%                 DATA(1,subject).DelayPupilData{j,1}=[forColumnAvg(j,:)];
%             end
            
                
        AllPupilDataOfInt=[AllPupilDataOfInt;(DATA(subject2).DelayPupilData)];
     
end


data.CleanedDelayData=AllPupilDataOfInt;

for i=1:height(data)
    
    data.MattMeanPupil{i,1}=nanmean(data.CleanedDelayData{i,1}(1,1000:3000));
    data.MattChangePupil{i,1}=nanmean(data.CleanedDelayData{i,1}(1,2000:3000))-nanmean(data.CleanedDelayData{i,1}(1,1000:2000));
   data.Baseline{i,1}=nanmean(data.CleanedDelayData{i,1}(1,900:1000));
    data.Change{i,1}=max(data.CleanedDelayData{i,1}(1,:))-data.Baseline{i,1};

end

end
% Look for systematic error in response based on Target Location:
%     figure(5)
%         scatter(data.TargAngle,data.ErrWRTOptimal);
%         title({'Systematic Thetal Error by Location'});
%         ylabel('Thetal Error (Degrees)');
%         xlabel('Targ Location (Degrees)');

        
        data.smoothedAvgRespErr=zeros(length(data.ErrWRTTarg),1);
        data.correctedErrWRTOptimal=zeros(length(data.ErrWRTTarg),1);
        smoothedAvgRespErr=smooth(data.TargAngle,data.ErrWRTOptimal, .5, 'rloess'); %May need to change span, comes before loess
        hold on
        [xx,ind] = sort(data.TargAngle);
%         plot(xx,smoothedAvgRespErr(ind),'r-')
        %RESIDUAL ERROR
        data.smoothedAvgRespErr=smoothedAvgRespErr;
        data.correctedErrWRTOptimal=data.ErrWRTOptimal-smoothedAvgRespErr;
%         legend('Trials','Loess fit');
%         axis([0,360,-40,40]);
       % ResponseVariablity=std(data{data{:,2}==1 & data{:,15}==2,24});
       
       
data=[data;ParseSubjectDataPeople( subject, pupil,plotit)];
       
       
       
end
