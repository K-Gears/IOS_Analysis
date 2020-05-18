function FitWhiskerStim(~)
%Fit whisker puff evoked hemodynamic response as sum of two gamma variate
%functions to capture dilation and overshoot seen in juvenile mice


StimStart=ProcessedData.dal_fr*5; %5 sec lead time for chunked data
StimEnd=StimStart+(5*ProcessedData.dal_fr);
TheData=ProcessedData.Contra_Puff.Run.Avg_Filt_Refl_barrels(StimStart:StimEnd);
Time=1:length(TheData);
FitParams=[-0.002,3,5,0.03,2,10,];%[-5,3,10,3,2,50];%Initial first guess of parameters
GammaFun=@(t,FitVars) (FitVars(1).*(t.^FitVars(2)).*exp(-t./FitVars(3)))+(FitVars(4).*(t.^FitVars(5)).*exp(-t./FitVars(6)));%Sum of to Gamma variate distribution functions to describe stimulus evoked dilation and overshoot
%GammaFun=@(t,FitVars) (FitVars(1).*(t.^FitVars(2)).*exp(-t./FitVars(3)));
ErrFun=@(t,FitVars)sum((GammaFun(t,FitVars)-TheData).^2);
options=optimset('PlotFcns',@optimplotfval);
options.MaxFunEvals=(1000*length(FitParams));
options.MaxIter=(1000*length(FitParams));
%[NewParams]=fminsearch(@(FitVars)ErrFun(Time,FitVars),FitParams,options);
A=[];
b=[];
Aeq=[];
beq=[];
lb=[-1,0,0,0,0,0];
ub=[0,Inf,Inf,1,Inf,Inf];
[NewParams]=fmincon(@(FitVars)ErrFun(Time,FitVars),FitParams,A,b,Aeq,beq,lb,ub,[],options);