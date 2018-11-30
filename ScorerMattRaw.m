function [ Scores ] = ScorerMattRaw( data )
%Gets Scores for a subject, must be done on parsed
delay=data.durationDelay;
type=data.sampleType;

        if type==1
    Scores=ceil(5+10*nanmean([data([data.correct]>.1).correct]));
        else
            Lowerbenchmark=mean(abs(diff([data.targetAngle])));
            Targets=[data.targetAngle];
            Means=[data.distMean];
            Means=[Means(2:end),Means(end)];
            Upperbenchmark=mean(abs(degAngDiff(Targets,Means)));
            
            PredictionError=mean(abs(degAngDiff([data.guessAngle],[data.NextAngle])));
            if PredictionError>Lowerbenchmark
                Payout=8;
            elseif PredictionError> 2/3*Lowerbenchmark+1/3*Upperbenchmark
                Payout=10;
            elseif PredictionError> 1/2*(Lowerbenchmark+Upperbenchmark)
                Payout=12;
            else
                Payout=15;
            end
            
           Scores=Payout;
       
        end
       
end

