function [A_agg,y_agg,t_agg,N_agg,p_refined_agg,r_agg] = aggregate_graph(A,y,t,p_refined,r,max_neigh_size)
%% This function aggregates the graph and all corresponding variables based
%% on the estimated cluster associations in the previous steps.
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
%
%Inputs:
%       A               : (numeric) An affinity/adjacency matrix of size
%                         N x N
%       y               : (numeric) A measurement vector of size N x 1
%       t               : (numeric) A time vector of size N x 1
%       p_refined       : (numeric) A refined partition vector of size
%                         N x 1. It includes cluster assignments that have
%                         been obtained via local moving of nodes and
%                         refinement steps in previous.
%       r               : (numeric) A trend vector of size N x 1. The
%                         vector has been updated based on local moving of
%                         nodes in the previous step.
%       max_neigh_size  : (numeric) The maximum number of neighbors whose
%                         edge weights that are attached to a considered
%                         cluster will be evaluated. The default setting
%                         considers all neighbors which is time consuming.
%                         A small value of the parameter is time efficient
%                         but it leads to small clusters in local moving
%                         of nodes step.
%
%Outputs :
%       A_agg           : (numeric) The aggregated affinity/adjacency
%                         matrix of size N_agg x N_agg
%       y_agg           : (numeric) The aggregated measurement vector of
%                         size N_agg x 1
%       t_agg           : (numeric) An aggregated time vector of size
%                         N_agg x 1
%       N_agg           : (numeric) Number of observations in the aggregated
%                         set. N_agg = number of clusters in p_refined
%       p_refined_agg   : (numeric) The aggregated refined partition vector
%                         of size N_agg x 1.
%       r_agg           : (numeric) The aggregated pattern vector of size
%                         N_agg x 1.
%
%
% version  : 24.05.2024
% Author   : Aylin Tastan
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

if nargin < 6 || isempty(max_neigh_size)
    max_neigh_size = length(p_refined)-1;
end

%Obtain the number of clusters and their indices
[~,ind_clusters,p_refined] = unique(p_refined,'stable');
N_agg = max(p_refined);

%Aggregate all variables
t_agg = t(ind_clusters);
r_agg = r(ind_clusters);
y_agg = y(ind_clusters);
p_refined_agg = p_refined(ind_clusters);

%Aggregate the affinity/adjacency matrix
ind_clusters = [ind_clusters;length(p_refined)+1];
A_agg = zeros(N_agg);

for i = 1:N_agg
    for j = i:N_agg
        if (j > i+max_neigh_size)
            break;
        else
            A_agg(i,j) = calc_matsum(A(ind_clusters(i):(ind_clusters(j+1)-1),ind_clusters(j):(ind_clusters(j+1)-1)));
        end
    end
end

A_agg = A_agg + tril(A_agg.',-1); %symmetrize

end