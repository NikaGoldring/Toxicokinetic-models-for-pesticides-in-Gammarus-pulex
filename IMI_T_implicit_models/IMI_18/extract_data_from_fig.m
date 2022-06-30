clear all;
close all;
clc;

%% Specify file for data extraction
open('20211213_IMI_standard_18.fig') ;
h = findobj(gcf, 'Type', 'line');

%% Extract xdata
xdata= get(h(4,1),'XData'); 

%% Extract ydata
%ydata= get(h,'YData');
ydata3= get(h(4,1),'YData'); % line of upper CI
ydata2= get(h(5,1),'YData'); % line of lower CI
ydata1= get(h(24,1),'YData'); % line of model fit

%% Export data to txt
fig01 = []; %create empty table
fig01(:,1) = xdata ; % xdata for all lines
fig01(:,2) = ydata3 ; % line of upper CI
fig01(:,3) = ydata2 ; % line of lower CI
fig01(:,4) = ydata1 ; % line of model fit
dlmwrite('20220125_IMI_standard_18.txt', fig01, ','); %write dataframe in txt
