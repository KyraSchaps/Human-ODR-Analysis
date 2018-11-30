function [ Errfits] = fMinLinearfit( data,MaxTACP, LargeCutOff, Cutoff, SessCut, endsession,Errsize,delays )

holdingData=data;
CUTOFF = LargeCutOff;% %45; %was 75
Lcp = [[holdingData.CPMagnitude]>CUTOFF];

Errfits(1:7,1)={'Intercept';'Delay'; 'TACPS'; 'Type'; 'Delay*Type'; 'Tacp*Type';'Delays*TACPS'};
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
           
for type=1:3
    
    
    if type<3
TrialsOfInterestANOVA= SessCondition & CPCondition &  holdingData.TACP<=MaxTACP & abs(holdingData.ErrWRTOptimal)<Errsize & holdingData{:,2}==type;
        yfit=@(x) x(1)+x(2)*DELAYS+x(3)*TACPS+x(4)*DELAYS.*TACPS;
    else
        TrialsOfInterestANOVA=  CPCondition &  holdingData.TACP<=MaxTACP & abs(holdingData.ErrWRTOptimal)<Errsize;
        yfit=@(x) x(1)+x(2)*DELAYS+x(3)*TACPS+x(4)*TYPE+x(5)*DELAYS.*TYPE+x(6)*TACPS.*TYPE+x(7)*DELAYS.*TACPS;

    end
y=abs(holdingData.ErrWRTOptimal(TrialsOfInterestANOVA));
MeanDiffs=holdingData.MeanTargDiff(TrialsOfInterestANOVA);
TACPS=holdingData.TACP(TrialsOfInterestANOVA);
TACPS(TACPS==0)=1;
DELAYS=holdingData.DelayTime(TrialsOfInterestANOVA);
TYPE=holdingData{TrialsOfInterestANOVA,2}-1;


% for i=1:length(delays)
% figure(20+type)
% subplot(1,2,i)
% hold on
% if type==1
%     response='Previous Sample Mean';
% elseif type==2
%     response='Generative Mean';
% end
% scatter(TACPS(DELAYS==delays(i)),y(DELAYS==delays(i)));
% % boxplot(y(DELAYS==delays(i)), TACPS(DELAYS==delays(i)));
% title(['Errors in Delay=' num2str(delays(i))]);
% ylabel(['Error from ',response]);
% xlabel(TACPS);
% end


% anovan(y,{TACPS,DELAYS},'model','interaction','varnames',{'TACPS','Delays'})


if type<3
            yfit=@(x) x(1)+x(2)*DELAYS+x(3)*TACPS+x(4)*DELAYS.*TACPS;

            
else
            yfit=@(x) x(1)+x(2)*DELAYS+x(3)*TACPS+x(4)*TYPE+x(5)*DELAYS.*TYPE+x(6)*TACPS.*TYPE+x(7)*DELAYS.*TACPS;

end

%Want to minimix (y real-yfit)^2  yfit=x*Means+(1-X)*Targs subject to
%constraint 0<x<1
fun = @(x)sum((y-yfit(x)).^2);
if type<3
x=fmincon(fun,[1,1,1,1],[],[],[],[],[],[]); 


else
x=fmincon(fun,[1,1,1,0,0,0,1],[],[],[],[],[],[]);
DELAYS=categorical(DELAYS);
TYPE=categorical(TYPE);
test=table(y,DELAYS,TACPS, TYPE, 'VariableNames',{'err','DELAYS','TACPS', 'TYPE'});
fit=fitlm(test,'err~DELAYS*TYPE*TACPS')
disp('done');
end
for i=1:length(x)
Errfits(i,type+1)={x(i)};
end

    
end

    Errfits=cell2table(Errfits,'VariableNames',{'Factor','Only_Per','Only_Pre','All'});

end
