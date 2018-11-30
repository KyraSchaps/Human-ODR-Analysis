function [ Scores ] = ScorerMatt( data,delays )
%Gets Scores for a subject, must be done on parsed



    for delay=1:length(delays)
        for type=1:2
        if type==1
    Scores(delay,type)=min(10,ceil(5+10*nanmean([data{[data{:,16}]>.02 & data.DelayTime==delays(delay) & data{:,2}==type,16}]')));
        else
            TrialsOfInt=[data{:,16}]>.02 & data.DelayTime==delays(delay) & data{:,2}==type;
            if sum(TrialsOfInt)>0
            Lowerbenchmark=mean(abs(diff(data{TrialsOfInt,3})));
            Targets=data{TrialsOfInt,3};
            Means=data{TrialsOfInt,5};
            Means=[Means(1,1);Means(1:end-1)];
            Upperbenchmark=mean(abs(degAngDiff(data{TrialsOfInt,3},Means)));
            PredictionError=mean(abs(degAngDiff(data{TrialsOfInt,7},data{TrialsOfInt,4})));
            disp(PredictionError);
            disp(Upperbenchmark);
            disp(Lowerbenchmark);
        
            if PredictionError>Lowerbenchmark
                Payout=5;
            elseif PredictionError> 2/3*Lowerbenchmark+1/3*Upperbenchmark
                Payout=8;
            elseif PredictionError> 1/2*(Lowerbenchmark+Upperbenchmark)
                Payout=9;
            else
                Payout=10;
            end
            
           Scores(delay,type)=Payout;
            else
                Scores(delay,type)=nan;
            end
        end
        end
    end

end

