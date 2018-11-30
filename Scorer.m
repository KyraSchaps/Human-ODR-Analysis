function [ Scores ] = Scorer( data )
%Gets Scores for a subject
num_subjects = size(unique(data{:,1}),1);
delays=[1,3,6];
for i=1:num_subjects
    for delay=1:3
        for type=1:2
    Scores(delay,type)=(1+.75*(type-1))*nanmean([data{[data{:,16}]>.1 & data.DelayTime==delays(delay) & data{:,2}==type,16}]');
        end
    end
end
end

