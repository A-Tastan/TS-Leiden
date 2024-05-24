function [p,r] = fast_and_local_moving_of_vertices(A,t,N,p,r,q,c_Tuk,T_weights)
%% This function performs Step 1 of TS-Leiden Algorithm.
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
% !NOTE :
%  This code requires an additional function 'madn', 'wtuk' and their all 
%  dependencies which are available in :
%
%  https://github.com/RobustSP/toolbox/tree/master/codes
%
%
%Inputs:
%       A         : (numeric) An affinity/adjacency matrix of size N x N
%       t         : (numeric) A time vector of size N x 1
%       N         : (numeric) Number of observations
%       p         : (numeric) A partition vector of size N x 1. It includes
%                   cluster assignments. If the vector is not defined, it
%                   assumes that all nodes are assigned to the different
%                   clusters.
%       r         : (numeric) A trend vector of size N x 1. If it is not
%                   specified, the algorithm assumes all trends are unknown
%                   which means that all variables are zero valued in this
%                   vector.
%       c_Tuk     : (numeric) Tuning constant for Tukey's weight function.
%                   The default value is 3.4437 for 85 % asymptotic relative
%                   efficiency.
%       T_weights : (numeric) Threshold for Tukey's weight function. It
%                   removes neighbors whose weight is smaller than the
%                   defined threshold. The default value is zero since
%                   Tukey's weight function gives zero weight to the
%                   determined outliers.
%
%Outputs:
%        p         : (numeric) The updated partition vector of size N x 1.
%                    Cluster assignments are obtained based on the local
%                    moving of nodes.
%        r         : (numeric) The updated trend vector of size N x 1. It
%                    assumes that the clusters share the same trend.
%
%
% version  : 24.05.2024
% Author   : Aylin Tastan
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

if nargin < 4 || isempty(p)
    p = 1:N;  %if the initial partition is not specified, it starts with a singleton partition
end

if nargin < 5 || isempty(r)
    r = zeros(N,1);  %if the initial trend vector is not specified, it starts with a zeros vector.
end

if nargin < 6 || isempty(q)
    q = 1:N;  %if nodes to visit vector is not specified, it includes all vertices
end

if nargin < 7 || isempty(c_Tuk)
    c_Tuk = 3.4437; %Tukey's tuning constant for real valued data
end

if nargin < 8 || isempty(T_weights)
    T_weights = 0;
end


%% Compute modularity and vertex to cluster edge weights
TOL = 1e-10;
%A = double(A); % convert adjacency matrix to double format

d = sum(A,"all");  % sum of edge weights for every node
B = (A-(sum(A,2)*sum(A,1))/d)/d; % modularity matrix of size N x N
if ~issymmetric(B)
    B = (B+B.')/2; % symmetrize modularity matrix
end

num_clusters = max(p);
H = zeros(num_clusters,num_clusters); % a matrix containing sum of vertex-to-module edge weights
for i = 1:max(p)  % loop over modules
    H(:,i) = sum(B(:,p==i),2);
end

%% Perform local moving for every vertex associated with an undetermined pattern
while ~isempty(q)

    %Determine the node to visit and associated parameters
    v_visit = q(1);  %node to visit
    q(1) = [];
    c_visit = p(v_visit); %cluster of the considered node
    p_visit = r(v_visit); %pattern of the considered node

    %Determine neighbors of visited node robustly
    [~,neigh_v_visit] = find(A(v_visit,:));
    neigh_v_visit(neigh_v_visit == v_visit) =[];
    weights = wtuk(abs(t(neigh_v_visit)-median(t(neigh_v_visit)))/madn(t(neigh_v_visit)),c_Tuk);
    neigh_v_visit= neigh_v_visit(weights > T_weights); %Define neighbors robustly

    %Calculate change in quality for neighbors of v_visit
    change_in_quality = H(v_visit,neigh_v_visit) - H(v_visit,c_visit) + B(v_visit,v_visit);

    %Determine the local move which provides the maximum quality improvement
    [max_quality_change,loc_max_quality_change] = max(change_in_quality); % maximal increase in quality and corresponding module

    %If the maximum quality improvement is positive, perform the local move
    if (max_quality_change > TOL)
        p(v_visit) = p(neigh_v_visit(loc_max_quality_change));
        r(v_visit) = r(neigh_v_visit(loc_max_quality_change));
    end
end


end
