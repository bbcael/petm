clear all; close all; clc; % clear workspace
c = readNPY('cenogrid_data.npy'); % load data
T = c(:,1); C = c(:,2); O = c(:,3); clear c; % label variables
w = 0.0302; % minimum time window size so PETM=1 time-point. 
% SENSITIVITY TEST #1 – adjust this parameter: up to 0.097 makes sense given data characteristics [n.b. doesn't change conclusions]
t = 0:w:-min(T); % define homogeneous time vector
i = find(t>55.8345); i = i(1); % find late edge of PETM time-point window
dw = 55.83445-t(i); t = -(t+dw); % find offset needed to keep PETM=1 time-point
for i = 1:(length(t)-1); % average d18O and d13C within time bins
    o(i) = nanmean(O(T<t(i) & T>t(i+1)-1e-12));
    % SENSITIVITY TEST #2 – use mean rather than median [n.b. doesn't change conclusions]
    c(i) = nanmean(C(T<t(i) & T>t(i+1)-1e-12));
end
t = (t(1:end-1)+t(2:end))./2; % make time-point middle of time-window
clear i dw C T O; % clear extra variables
os = smoothdata(o,'movmedian',8); cs = smoothdata(c,'movmedian',8); % define moving mean – window size identified with parcorr
% SENSITIVITY TEST #3 – use moving mean rather than median [n.b. doesn't change conclusions]
% SENSITIVITY TEST #4 – increase smoothing window [n.b. doesn't change conclusions]
oa = o-os; ca = c-cs; oa(oa==0) = NaN; ca(ca==0) = NaN; % define fluctuations
xo = sort((oa(oa>-max(oa)))); xc = sort((ca(ca>-max(ca)))); % define data for fitting distributions – may be one- or two-sided
% SENSITIVITY TEST #5 – try fit for e.g. different ranges of "os" for state-dependence of exponential parameter [n.b. weak state-dependence, but doesn't change conclusions]
[eo,uo] = mle(xo,'pdf',@(x,b)1/(2*b)*exp(-abs(x)/b),'Start',[.1]); % empirical estimate of laplace parameter 
% SENSITIVITY TEST #6 – try different probability distributions [n.b. Laplace is best fit with fewest free parameters and simplest phenomenological explanation, doesn't change conclusions]
ceo = heaviside(xo)-.5.*sign(xo).*exp(-abs(xo)/eo(1)); y = linspace(0,1,length(xo)); vo = max(ceo-y) + max(y-ceo); % find empirical-theoretical distribution mismatch
[ec,uc] = mle(xc,'pdf',@(x,b)1/(2*b)*exp(-abs(x)/b),'Start',[.1]); % same as above for d13C
cec = heaviside(xc)-.5.*sign(xc).*exp(-abs(xc)/ec(1)); y = linspace(0,1,length(xc)); vc = max(cec-y) + max(y-cec);  % same as above for d13C
clear cec ceo xc xo y; % clear extra variables
uc = (uc(2,:)-uc(1,:))./3.92; uo = (uo(2,:)-uo(1,:))./3.92; % convert to 1-sigma
nb = 10000; % bootstrap iterations
for i = 1:nb; % generate probability distribution of maxima for Cenozoic
    xo = exprnd(eo+randn(1).*uo,1,round(6.6e7/(w.*1e6))); % simulate fluctuations
    % SENSITIVITY TEST #7 -- increase 6.6e7 (duration of Cenozoic) to age of Earth (4.5e9)
    xc = exprnd(ec+randn(1).*uc,1,round(6.6e7/(w.*1e6))); % simulate fluctuations
    Omax(i) = max(xo); Cmax(i) = max(xc); clear xo xc; % find largest fluctuations
end
clear i w; % clear extra variables
xo = exprnd(eo,nb,length(oa)); % generate fake data to test if sample size can explain distribution fit error
for i = 1:nb; % fit generated data with generating distribution
    [y,x] = ecdf(xo(i,:)); yy = 1-exp(-x./eo); vob(i) = max(y-yy) + max(yy-y); % find theoretical-sampled distribution mismatch
end
xc = exprnd(ec,nb,length(ca)); % same for d13C
for i = 1:nb;
    [y,x] = ecdf(xc(i,:)); yy = 1-exp(-x./ec); vcb(i) = max(y-yy) + max(yy-y); % same for d13C
end
% output some results:
sum(Omax>-min(oa))./length(Omax) % probability of a PETM-sized d18O fluctuation during the Cenozoic
sum(Cmax>-min(ca))./length(Cmax) % probability of a PETM-sized d13C fluctuation during the Cenozoic
sum(vob>vo)./length(vob) % probability of Laplace-distributed data having as much of a sample-distribution mismatch as the d18O observations 
sum(vcb>vc)./length(vcb) % probability of Laplace-distributed data having as much of a sample-distribution mismatch as the d13C observations 
clear ans nb xo i y x yy xc; % clear extra variables

%% theory: d(delta')/dt = f - r*sign(delta')

clc; % clear command line

x(1,1:length(c)) = 0; % dummy variable -- delta prime in paper notation
f = randn(1000,length(c)); % gaussian forcing

for i = 1:size(f,1); % simulate random forcing - constant relaxation dynamics
    x(i+1,:) = x(i,:) + .01.*f(i,:) - .0007.*sign(x(i,:));
end % do 10,000 iterations of equal length to observations

x = x(2:end,:); x = x(:); % remove initial timestep
b = mle(x,'pdf',@(x,b)1/(2*b)*exp(-abs(x)/b),'Start',[.1]); % estimate laplace parameter

[y,xx] = ecdf(x); yy(xx<0) = .5*exp(xx(xx<0)./b); yy(xx>0) = 1-.5.*exp(-xx(xx>0)./b); v = max(abs(y-yy')) + max(abs(yy'-y)); % find model-laplace mismatch
[v b] % report Kuiper statistic and b parameter value of fit

clear f i b y xx yy v; % clear extra variables

%% figure 1a

close all; clc;
figure
subplot(121)
[hoy,hox] = histcounts(oa(oa>min(oa)),-.495:.03:.5, 'Normalization', 'probability'); hox = (hox(1:end-1)+hox(2:end))./2;
[hcy,hcx] = histcounts(ca(ca>min(ca)),-.495:.03:.5, 'Normalization', 'probability'); hcx = (hcx(1:end-1)+hcx(2:end))./2;
[hty,htx] = histcounts(x(x>min(x)),-.495:.03:.5, 'Normalization', 'probability'); htx = (htx(1:end-1)+htx(2:end))./2;
p4 = scatter(htx,hty,200,[.6 .05 .4],'s','filled');
hold on;
p3 = scatter(hcx,hcy,200,[0.05 .5 .5],'filled');
p2 = scatter(hox,hoy,200,[0.7 0.4 0.05],'^','filled');
x = -.5:.001:.5;
l = .0708.*exp(-abs(x)./.0708); l = 29.*l./sum(l);
p1 = plot(x,l,'k','linewidth',2);
clear hoy hox hcx hcy x;
axis([-.4 .4 0 .22])
box on;
lgnd = legend([p1 p2 p3 p4],'Laplace PDF','\delta^{18}O Anomalies','\delta^{13}C Anomalies','CRR Model');
set(lgnd,'fontsize',15,'location','northwest')
set(gca,'fontsize',15,'ytick',[],'xtick',-.4:.1:.4)
xlabel(['Anomaly (', char(8240),')'])
ylabel('Probability Density')
clear p1 p2 p3 lgnd;
text(.345,max(l),'A','fontsize',30)
clear l;

%% figure 1b

subplot(122)
[hoy,hox] = histcounts(Omax,min(Omax):.01:max(Omax), 'Normalization', 'probability'); hox = (hox(1:end-1)+hox(2:end))./2;
[hcy,hcx] = histcounts(Cmax,min(Omax):.01:max(Omax), 'Normalization', 'probability'); hcx = (hcx(1:end-1)+hcx(2:end))./2;
hcy = smoothdata(hcy,'gaussian',15); hoy = smoothdata(hoy,'gaussian',15);
p2 = plot(hcx,hcy,'linewidth',2,'color',[0.05 .5 .5]);
hold on;
p1 = plot(hox,hoy,'color',[0.7 0.4 0.05],'linewidth',2);
clear hox hcx hcy;
axis([min(Omax-.05) 2.5 0 .05])
box on;
set(gca,'fontsize',15,'ytick',[]);
xlabel(['Largest Anomaly in 66 Ma (', char(8240),')'])
p3 = plot(linspace(max(abs(ca)),max(abs(ca))),linspace(-1,1),'--','linewidth',2,'color',[.05 .5 .5])
p4 = plot(linspace(max(abs(oa)),max(abs(oa))),linspace(-1,1),'--','linewidth',2,'color',[.75 .4 .05])
[~,ic] = max(abs(ca)); ca((ic-1):ic) = []; oa(ic) = []; t(ic) = []; % remove PETM to find next-largest anomalies
oa(t>-4) = []; % remove G-IG d18O anomalies because these are qualitatively different dynamics
% SENSITIVITY TEST #8 -- don't remove G-IG d18O anomalies in to see that even these aren't larger than expected
% SENSITIVITY TEST #9 -- remove G-IG cycle (last 4 Ma) entirely from analysis [n.b. does not affect conclusions]
p5 = plot(linspace(max(abs(ca)),max(abs(ca))),linspace(-1,1),':','linewidth',2,'color',[.05 .5 .5])
p6 = plot(linspace(max(abs(oa)),max(abs(oa))),linspace(-1,1),':','linewidth',2,'color',[.75 .4 .05])
lgnd = legend([p1 p2 p4 p3 p6 p5],'\delta^{18}O Expected Maximum Distribution','\delta^{13}C Expected Maximum Distribution','\delta^{18}O PETM Signature','\delta^{13}C PETM Signature','\delta^{18}O Next-Largest Anomaly','\delta^{13}C Next-Largest Anomaly');
set(lgnd,'fontsize',15,'location','northeast')
text(min(Omax-.01),.97.*max(hoy),'B','fontsize',30)
clear hoy;