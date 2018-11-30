function [ fakechoices,Probs] = GenFakeChoices(LLR,H,noise,lapse,ChoiceNoiseFlag )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
% 
%0 means choose T2, 1 means chose T1


randoms=rand(length(LLR),1);
[Probs,Belief]=PredictedProb(LLR,H,noise,lapse );

if ChoiceNoiseFlag
fakechoices(randoms<Probs)=1;
fakechoices(randoms>=Probs)=0;
fakechoices=fakechoices';
else
    fakechoices=sign(Belief);
    fakechoices(fakechoices<0)=0;
end

end

