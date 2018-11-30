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
title(['type', num2str(j-1), 'CP cutoff ', num2str(cutoff),' ' , ' Trial Cutoff ', num2str(endsession)]);
hold on


for i=1:length(subjects)
hold on
Baseline=[Baseline;cell2mat(data.Linfit_Pup_Mean{i,1}{:,j})];
subjectnames=[subjectnames;repmat(data.SubjectFile(i),4,1)];
item=[item;[1:4]'];
end

gscatter(item,Baseline,subjectnames);
end

