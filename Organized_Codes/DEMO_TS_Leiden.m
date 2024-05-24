%% This demo file runs TS-Leiden Algorithm for risky time determination of
%% hypoglycemia. For details, see :
%
%   A. Tastan, C. Escorihuela-Altaba, J. Garcia-Tirado and K. Riesen,
%  "Clustering Time Series Data for Personalized Type 1 Diabetes
%  Management", in Proc. Int. Conf. Pattern Recognit. Artif. Intell.
%  (ICPRAI2024), 2024.
% 
%  Copyright (C) 2024 Aylin Tastan. All rights reserved.
% 
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with this program. If not, see <https://www.gnu.org/licenses/>.
%
%
% !NOTE :
% This code requires an additional function 'madn', 'wtuk' and their all 
% dependencies which are available in :
%
% https://github.com/RobustSP/toolbox/tree/master/codes
%
% This code requires an additional function 'fast_NVG' and its all 
% dependencieswhich is available in :
%
% https://ch.mathworks.com/matlabcentral/fileexchange/70432-fast-natural-visibility-graph-nvg-for-matlab
%
% version : 24.05.2024
% author : Aylin Tastan
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
warning off
plotting = 0;

%% Call CGM recording and format as time series
Data =  table2cell(readtable('cgm.csv'));
datetime_cgm = Data(:,2);
for i = 1:(length(datetime_cgm)-1)
   diff_cgm(i) = minutes(datetime_cgm{i+1,1}-datetime_cgm{i,1});
end
time_cgm = [minute(datetime_cgm{1,1}),minute(datetime_cgm{1,1}) + cumsum(diff_cgm)].';
observation_cgm = cell2mat(Data(:,3));


%% Perform the TS-Leiden Algorithm
tic
% for Mc = 1:100
[A_cgm,y,t,p,r] = perform_TS_Leiden_Algorithm(observation_cgm,time_cgm);
% end
toc

%Plot trend determination results
if plotting
    plot_trend_determination_results(observation_cgm,time_cgm,t,[r;0]);
end

%Summarize TS-Leiden Algorithm-based clustering and trend determination results
p_all = zeros(length(time_cgm),1);
t =[t;time_cgm(end)+1];
for i = 1:length(p)
       ind_clust_i = find(time_cgm>=t(i) & time_cgm <t(i+1));
       p_all (ind_clust_i) = i;
       r_all (ind_clust_i) = r(i);
end

results = [observation_cgm,p_all,r_all.'];

%% Risky time determination for hypoglycemic events
%Find clusters associated with decreasing trends that lead to hypoglysemia
[clusters,ind_endpoint,~] = unique(results(:,2),'last','legacy');
datetime_hypoglysemia = [];
for j = 1:length(clusters)
    if(j == 1)
        if(results(ind_endpoint(j),3) == -1 & results(ind_endpoint(j),1) < 70)
            datetime_hypo = [datetime_cgm{1,1},datetime_cgm{ind_endpoint(j),1}];
            datetime_hypoglysemia = [datetime_hypoglysemia;datetime_hypo];
        end
    else
        if(results(ind_endpoint(j),3) == -1 & results(ind_endpoint(j),1) < 70)
            datetime_hypo = [datetime_cgm{ind_endpoint(j-1)+1,1},datetime_cgm{ind_endpoint(j),1}];
            datetime_hypoglysemia = [datetime_hypoglysemia;datetime_hypo];
        end
    end
end

%Observe time information for hypoglysemia
time_hypoglysemia = [];
for i = 1:length(datetime_hypoglysemia)
    time_hypo = [timeofday(datetime_hypoglysemia(i,1)),timeofday(datetime_hypoglysemia(i,2))];
    time_hypoglysemia = [time_hypoglysemia;time_hypo];
end

[time_hypoglysemia(:,1),ind_hypo] = sort(time_hypoglysemia(:,1));
time_hypoglysemia(:,2) = time_hypoglysemia(ind_hypo,2);

%Segment hypoglysemia times based on sin and cos features
min_timehypo = minutes(time_hypoglysemia(:,1));
for j = 1:length(time_hypoglysemia)
    hypoglysemia_fea(j,1) = sin((min_timehypo(j)*2*pi)/1440);
    hypoglysemia_fea(j,2) = cos((min_timehypo(j)*2*pi)/1440);
end

c_hat = kmedoids(hypoglysemia_fea,3,'Replicates',10);

time_hypoglysemia_ord = [];
for k = 1:max(c_hat)
    idx_k = find(c_hat==k);
    time_label_pair = [time_hypoglysemia(idx_k,:),k*ones(length(idx_k),1)];
    time_hypoglysemia_ord = [time_hypoglysemia_ord;time_label_pair];
end

















