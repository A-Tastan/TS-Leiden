function [p_refined,r,num_outliers] = refinement_of_the_partition(y,p,r,min_num_samples)
%% This function performs Step 2 of TS-Leiden Algorithm.
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
%       y              : (numeric) A measurement vector of size N x 1
%       p              : (numeric) A partition vector of size N x 1. It
%                        includes cluster assignments that have been
%                        obtained via local moving of nodes in the previous
%                        step.
%       r              : (numeric) A trend vector of size N x 1. The
%                        vector has been updated based on local moving of
%                        nodes in the previous step.
%      min_num_samples : (numeric) The minimum number of samples that is
%                        necessary to form a cluster. The default value is
%                        one based on singleton clusters.
%
%Outputs :
%        p_refined     : (numeric) The refined partition vector of size
%                        N x 1.
%        r             : (numeric) The updated trend vector of size N x 1.
%        num_outliers  : (numeric) Number of outliers estimate
%
%
% version  : 24.05.2024
% Author   : Aylin Tastan
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

if nargin < 4 || isempty(min_num_samples)
    min_num_samples = 1;
end

%Determine number of clusters and order cluster assignments
[p_unique,~,p_refined] = unique(p,'stable');

%Split node subsets based on their trends
[p_refined,r,ind_outliers,num_clusters] = split_vertex_subsets(y,p_refined,r,length(p_unique));

%Refine trends correspond to outliers
[p_refined,r,num_outliers] = refine_outliers(y,p_refined,r,ind_outliers,num_clusters);

end%endfunction

