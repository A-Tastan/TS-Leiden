function [p_refined,r,ind_outliers,num_clusters] = split_vertex_subsets(y,p_refined,r,num_clusters)
%% This function split clusters into different subsets based on their hidden
%% trends.
%
%  For details, see:
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
%Inputs:
%       y               : (numeric) A measurement vector of size N x 1
%       p_refined       : (numeric) A partition vector of size N x 1 which
%                         will be refined at this step. It includes cluster
%                         assignments that have been obtained via local
%                        moving of nodes in the previous step.
%       r               : (numeric) A trend vector of size N x 1. The
%                         vector has been updated based on local moving of
%                         nodes in the previous step.
%       num_clusters    : (numeric) Number of clusters in p_refined
%       min_num_samples : (numeric) The minimum number of samples that is
%                        necessary to form a cluster. The default value is
%                        one based on singleton clusters.
%
%Outputs :
%        p_refined      : (numeric) The refined partition vector of size
%                         N x 1.
%        r              : (numeric) The updated trend vector of size N x 1.
%        ind_outliers   : (numeric) A vector containing indices of
%                         estimated outliers. It's size is num_outliers x 1
%        num_outliers   : (numeric) Number of outliers estimate
%
%
% version  : 24.05.2024
% Author   : Aylin Tastan
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

ind_outliers =[];

%% Split cluster into sequential trends and determine outliers
for i = 1:p_refined(end)
    ind_cluster_i = find(p_refined==i); %indice of nodes that are associated with cluster i

    %Split cluster based on the hidden trends
    if (length(ind_cluster_i) > 1) %if not a singleton cluster
        trend_vec = sign(diff(y(ind_cluster_i)));
        trend_vec = [trend_vec;trend_vec(end)];


        %Change trend information of nodes that are correspond to a positive trend
        pos_trend = ind_cluster_i(find(trend_vec >= 0));
        if ~isempty(pos_trend)
            pos_size = length(pos_trend);
            pos_sequential = pos_trend(strfind(([pos_size;diff(pos_trend);pos_size]==1).',[0 1]):strfind(([pos_size;diff(pos_trend);pos_size]==1).',[1 0]));
            if ~isempty(pos_sequential)
                r(pos_sequential) = 1 ;
            end
        else
            pos_sequential = [];
        end

        %Assign nodes correspond to a negative trend into a different cluster and update the trend vector
        neg_trend = ind_cluster_i(find(trend_vec < 0));
        if ~isempty(neg_trend)
            neg_size = length(neg_trend);
            neg_sequential = neg_trend(strfind(([neg_size;diff(neg_trend);neg_size]==1).',[0 1]):strfind(([neg_size;diff(neg_trend);neg_size]==1).',[1 0]));
            if ~isempty(neg_sequential)
                num_clusters = num_clusters + 1;
                p_refined(neg_sequential) = num_clusters;
                r(neg_sequential) = -1 ;
            end
        else
            neg_sequential =[];
        end


        ind_outliers = [ind_outliers;setdiff(ind_cluster_i,[pos_sequential;neg_sequential])];

    end
end

end%endfunction

