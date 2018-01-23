%% Restitution generator %%
% Description: Uses the oapRestitution.m function to calculate all action
% potential values and assemble a restitution curve. It can be modified to
% generate a dose response curve as well.
%
% Input:
%   None
%
% Output:
%   None
%
% Author: Christopher Gloschat
% Date: August 31, 2016
%
%
% Modification Log:
%
%
%% Code %%
% cycle length or drug concentration
% cycle = [300 275 250 225 200 175];
cycle = [0 1 2 5 10 20 40 60];
region = 'RV'; 
% preallocste variables
x = [];
y = [];
std = [];
num = [];
% save file
% tmp = sprintf('save %s_Restitution x y std num',region);
tmp = sprintf('save %s_doseResponse x y std num',region);
eval(tmp)
for n = 1:length(cycle)
   % find average APD
   oapRestitution(cycle(n),region)
   % load results
   load('Transfer.mat')
   % load restitution variables
   tmp = sprintf('load(''%s_doseResponse.mat'')',region);
   eval(tmp)
   % save out transfer variables
   x = [x conc];
   y = [y apd80TimeAve];
   std = [std apd80TimeStd];
   num = [num numPeaks];
   % clear transfer variables
   clear conc apd80TimeAve apd80TimeStd numPeaks
   % save restitution values
   tmp = sprintf('save %s_doseResponse x y std num',region);
   eval(tmp)
end

% visualize
% % % figure
% % % errorbar(x,y,std,'ko--')
% % % tmp = sprintf('title(''%s Dose Response @ 300ms Pacing'',''FontSize'',22)',region);
% % % eval(tmp)
% % % xlabel('Concentration (uM)','FontSize',20)
% % % ylabel('APD_8_0 (ms)','FontSize',20)
% % % set(gca,'FontSize',18)
% visualize
IC50 = log10(40e-6); % IC50 of 30 um
xb = -5.8:0.001:-4.2;
yb = repmat(min(y)+(max(y)-min(y)),[size(xb,1) size(xb,2)])./(1+exp((IC50-xb)/-0.35));
figure
errorbar(log10(x(3:end)*1e-6),y(3:end),std(3:end),'ko')
hold on
plot(xb,yb,'k--')
[~,IC50y] = min(abs(xb-IC50));
IC50y = yb(IC50y);
plot([IC50 IC50],[0 IC50y],'k')
plot([-5.8 IC50],[IC50y IC50y],'k')
hold off
tmp = sprintf('title(''%s Dose Response @ 300ms Pacing'',''FontSize'',22)',region);
eval(tmp)
xlabel('Concentration (uM)','FontSize',20)
ylabel('APD_8_0 (ms)','FontSize',20)
set(gca,'XLim',[-5.8 -4.2],'YLim',[20 180])
set(gca,'FontSize',18)


