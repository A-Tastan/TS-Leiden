function [p_refined,r,num_outliers] = refine_outliers(y,p_refined,r,ind_outliers,num_clusters)
%% This function recovers undesired effects of outliers based on their actual
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
%       ind_outliers    : (numeric) A vector containing indices of
%                         estimated outliers. It's size is num_outliers x 1
%       num_clusters    : (numeric) Number of clusters in p_refined
%
%Outputs :
%        p_refined      : (numeric) The refined partition vector of size
%                         N x 1.
%        r              : (numeric) The updated trend vector of size N x 1.
%        num_outliers   : (numeric) Number of outliers estimate
%
%
% version  : 24.05.2024
% Author   : Aylin Tastan
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


% Assign outliers into different singleton clusters and determine their trends
act_trend = sign(diff([y;y(end)]));
ind_outliers = sort(ind_outliers);
len_p_refined = length(p_refined);
num_outliers = length(ind_outliers);
if ~isempty(ind_outliers)

    %Define trend values correspond to outliers
    r(ind_outliers) = act_trend(ind_outliers); %determine their trends based on their neighbors;

    %Assign first outlier to a cluster
    if(ind_outliers(1)==1)
        if(abs(r(1)-r(2))<=1)

            p_refined(1) = p_refined(2);
        else
            num_clusters = num_clusters +1;
            p_refined(1) = num_clusters;
        end
    end

    %Assign last outlier to a cluster
    if(ind_outliers(end) == len_p_refined)
        if(abs(r(len_p_refined)-r(len_p_refined-1))<=1)

            p_refined(len_p_refined) = p_refined(len_p_refined-1);
        else
            num_clusters = num_clusters +1;
            p_refined(len_p_refined) = num_clusters;
        end
    end

    %Assign remaining outlier indices
    for i = 2:(length(ind_outliers)-1)

        if(abs(r(ind_outliers(i)) - r(ind_outliers(i)-1))<=1)

            p_refined(ind_outliers(i)) = p_refined(ind_outliers(i)-1);

        else if(abs(r(ind_outliers(i)) - r(ind_outliers(i)+1))<=1)

                p_refined(ind_outliers(i)) = p_refined(ind_outliers(i)+1);

        else
            num_clusters = num_clusters +1;
            p_refined(ind_outliers(i)) = num_clusters;

        end%endelseif

        end%endif

    end

end


end%endfunction