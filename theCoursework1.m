%% Task 1&2
clc
clear all
%% import data 
load('series1');
%% delete NaN observations
miss = isnan(logret);
logret = logret(~miss);  
dat = dat(~miss,:);
%%
date_num=datenum(dat);
%% calculate statistics
[mean,var,std,skew,kurt] = SummaryStats(logret);
%% implement JB-test
h1 = jbtest(logret);
% series1 is from a normal distribution
%% plot the time series  against timeline horizon.
nobs = size(logret,1);
seq = (1:nobs)';
year = linspace(1996.1, 2014.12,nobs);
plot(year(seq), logret);
xlim([1996,2014]);
title('logret from 1996 to 2014');
%% implement ADF test
h2 = adftest(logret);
% series1 is a stationary (no unit root)
%% implement ACF PACF
autocorr(logret);
% autocorrelation cuts off after 3 lags
parcorr(logret);
% partial autocorrelation dies off slowly
% series1 is a MA(3) process

%% Identify the best lag order q using BIC and AIC
%  ARMA(p,q) model with p,q = 0,1,...,maxlag

logret=logret./100; 
maxlag     = 4;    % We search over lag orders Q = 0,1,...,4
save_AIC   = zeros(maxlag+1,1);    % Save AIC for all possible lags of Q
save_BIC   = zeros(maxlag+1,1);    % Save BIC for all possible lags of Q
T=length(logret);
  
    for i = 0:maxlag
        Q=(1:i);
 

            [b, LL, errors, SEregression, di] = armaxfilter(logret,1,[],Q,[],[],[],maxlag);
            % we set up P=0, so as to test AIC, BIC for different MA(Q)
            % models
            save_AIC(i+1,1) = di.AIC;
            save_BIC(i+1,1) = di.SBIC;
       
    end
 %% conduct rolling window forecasting procedure on the MA(4) model  
 %identify the location of the first obs in every month
ind_0=zeros(T,1);
obs_0=zeros(T,1);
for i=2:T;          
    ind_0(i)=dat(i,2)-dat(i-1,2); % check if observations belong to same month.
    if ind_0(i)~=0; 
       obs_0(i)=length(ind_0(1:i,1));
    end 
end 
obs_0=obs_0(obs_0~=0); % delete zeros terms
obs_0=[1;obs_0;T]; % obs_1: indicate the location for the first observation in every month;

obs_n=obs_0(2:end)-obs_0(1:end-1); % calculate the obs number in each month and save as obs_n
obs_n=obs_n(obs_n~=0);  % delete zeros terms from obs_n

%% Conditional Forecasting with Rolling Estimation Window
% Ist step: determine the initial estimation window
% Initial estimation window: 2000/1/3-2015/12/31
T0=datenum('2008/1/31','yyyy/mm/dd');
T1=find(date_num==T0); % Locate the ending date for estimation window
T2=T1+1; % locate the starting date for forecasting window
T3=length(obs_0)-find(obs_0==T2); %  months in forecast window.


% for rolling estimation window, save updated parameters' estimations in 
% parameter_save

P=0;
Q=4;
par=zeros(P+Q+1,T3);

yhat_save=[];
mse=zeros(T3,1);
mae=zeros(T3,1);
temp=find(obs_0==T2); 
for i=1:T3; % loop for estimation window (in a monthly frequency)
       estwin_sta=obs_0(i); % estimation window starts
       estwin_end=obs_0(temp)-1; % estimation window ends
       est=logret(estwin_sta:estwin_end); % set up estimation window
       est=est-mean(est); % return- mean 
       [parameters,LL,errors]=armaxfilter(est,1,(1:P),(1:Q));

       
       par(:,i)=parameters; % save estimated parameters
       al=par(1,i);
       ph=par(2:Q+1,i);
       th=par(Q+2:end,i);
       forday=obs_m(temp); % number of forecast days in a month.
       % generate forecasting for next month.
       yf = armaforecast(est,forday,P,Q,al,ph,th);
       yhat_save=[yhat_save;yf(length(est)+1:end)];
       mse(i,1)=sum((logret(estwin_end+1: estwin_end+forday)-yf(length(est)+1:end)).^2)/forday;
       mae(i,1)=sum(abs(logret(estwin_end+1: estwin_end+forday)-yf(length(est)+1:end)))/forday;
       temp=temp+1;
end
%% Using Ljung-Box test
L=20; % Set up 20 lags for calculating the serial correlations among error terms.
Q_LB=ljungbox(logret,L) % Q=boxljung(x,h) gives Box Ljung statistics


%% Task 3
clc
clear all

%-----------import data--------------------------------%
load('DJIA.mat');  % import data

%% delete NaN observation
miss=isnan(price);
price=price(~miss);  
date=date(~miss);
%% calculate the log-returns for asset prices
% the log-return can be calculated as r_t=ln(p_t)-ln(p_(t-1))
% This way r_1=0
ret=[0;diff(log(price))].*100;

%% plot the curve for both return and price level
figure;
subplot(2,1,1);
nobs=size(ret,1);
seq=(1:nobs)';
year=linspace(2011.1, 2019.12,nobs);
plot(year(seq), ret);
xlim([2011,2019]);

subplot(2,1,2);
plot(year(seq), price);
xlim([2011,2019]);
title('Return from 2011 to 2019');

% represent dates of data_return as vector format;
datevector=datevec(date); 
% convert date into a numerical format;
date_num=datenum(date,'yyyy/mm/dd');
% length of return
T=length(ret);

%% identify the location of the first obs in every month
ind=zeros(T,1);
obs_1=zeros(T,1);
for i=2:T;          
ind(i)=datevector(i,2)-datevector(i-1,2); % check if observations belong to same month.
    if ind(i)~=0; 
       obs_1(i)=length(ind(1:i,1));
    end 
end 
obs_1=obs_1(obs_1~=0); % delete zeros terms
obs_1=[1;obs_1;T]; % obs_1: indicate the location for the first observation in every month;

obs_m=obs_1(2:end)-obs_1(1:end-1); % calculate the obs number in each month and save as obs_m
obs_m=obs_m(obs_m~=0);  % delete zeros terms from obs_m

%IND=cumsum(ind)+1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% estimate the full sample for ht
eps=ret-mean(ret);
[~, ~, ht_full] = tarch(eps,1,0,1,'NORMAL' );

%%  Conditional volatility forecasting with rolling window GARCH(1,1)
% 1st step: determine the initial estimation window
% initial estimation window: 2011/1/3-2019/12/30
T0=datenum('2015/12/31','yyyy/mm/dd');
T1=find(date_num==T0); % Locate the ending date for estimation window
T2=T1+1; % locate the starting date for estimation window
T3=length(obs_1)-find(obs_1==T2); %  months in forecast window.

% for each loop, save estimated ht as into  matrix ht_save;
%ht_save=zeros(T,T3);

% for rolling estimation window, save updated parameters' estimations in 
% parameter_save
par=zeros(3,T3);

mse=zeros(T3,1);
mae=zeros(T3,1);

% then for each loop, save new generated forecasts into matrix hf_save;
hf_save=zeros(1,1);
temp=find(obs_1==T2); 
for i=1:T3; % loop for estimation window (in a monthly frequency)
       estwin_sta=obs_1(i); % estimation window starts
       estwin_end=obs_1(temp)-1; % estimation window ends
       est=ret(estwin_sta:estwin_end); % set up estimation window
       est=est-mean(est); % return- mean 
       [parameters, LL, ht] = tarch(est,1,0,1,'NORMAL' );
       
       par(:,i)=parameters; % save estimated parameters
       ht_save(1:length(ht),i)=ht; % save estimated ht
       forday=obs_m(temp); % number of forecast days in a month.
       % generate forecasting for next month.
       ht_f=Garchforecast(est,ht,forday,parameters);
       
       hf_save=[hf_save;ht_f];
       hf_save=hf_save(hf_save~=0);
       mse(i,1)=sum((ht_full(estwin_end+1: estwin_end+forday)-ht_f).^2)/forday;
       mae(i,1)=sum(abs(ht_full(estwin_end+1: estwin_end+forday)-ht_f))/forday;
       temp=temp+1;
       
end

ht_compare=ht_full(T2:end-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare the forecasting results of hf with restimated ht
figure;
nobs=size(hf_save,1);
seq=(1:nobs)';
year=linspace(2016.1, 2019.12,nobs);
plot(year(seq), hf_save);
xlim([2016,2019]);
hold on;
plot(year(seq), ht_compare);
legend('Forecasting','Estimation')

%% identify the location of the first obs in every quarter
ind_quar=zeros(T,1);
obs_quar=zeros(T,1);
for i=2:T;          
ind_quar(i)=datevector(i,2)-datevector(i-1,2); % check if observations belong to same quarter.
    if ind_quar(i)~=0 & datevector(i-1,2)==3 | datevector(i-1,2)==6 | datevector(i-1,2)==9 | datevector(i-1,2)==12;
       obs_quar(i)=length(ind(1:i,1));
    end 
end 
