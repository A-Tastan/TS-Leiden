function [A,y,t,p,r] = perform_TS_Leiden_Algorithm(y,t)
%% This function performs TS-Leiden Algorithm. 
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
%       y        : (numeric) A measurement vector of size N x 1
%       t        : (numeric) A time vector of size N x 1
%
%Outputs:
%       A        : (numeric) The aggregated affinity/adjacency matrix of
%                   size N_agg x N_agg (N_agg = number of clusters in p)
%       y        : (numeric)The aggregated measurement vector of size
%                   N_agg x 1
%       t        : (numeric)The aggregated time vector of size N_agg x 1
%
%       p        : (numeric)The partition vector estimate of size N_agg x 1
%
%       r        : (numeric) The trend vector estimate of size N_agg x 1
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
% version  : 24.05.2024
% Author   : Aylin Tastan
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%% Step 0: Initialization
N = length(y);
N_out = N;
TOL = 5;
p = (1:N).'; %singleton clusters
r = zeros(N,1); %trends assumed to be unknown
q = 1:N;  %all vertices to be visited
c_Tuk = 3.4437; %Tukey's tuning constant for real valued data
T_weights = 0;
max_neigh_size = 10;
A = full(fast_NVG(y,t,'u',0)); %NVG


%% Repeating Steps (1-3) of TS-Leiden

while(1)

    %% Step 1: Fast and Robust Local Moving of Vertices
    [p,r] = fast_and_local_moving_of_vertices(A,t,N,p,r,q);

    %% Step 2: Refinement of the Partition
    [p,r,~] = refinement_of_the_partition(y,p,r);

    %% Step 3 : Aggregation of the Graph based on the Refined Partition
    [A,y,t,N,p,r] = aggregate_graph(A,y,t,p,r,max_neigh_size);

    q = find(r==0);

    if ~(length(q) > TOL)
        break;
    end

end

end