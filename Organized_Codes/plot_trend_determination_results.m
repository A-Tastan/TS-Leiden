function plot_trend_determination_results(y,x,t,r)
%% This function plots trend determination results.
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
% version : 24.05.2024
% author : Aylin Tastan
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
ind_pos = [];
ind_neg = [];
ind_out = [];
for i = 1:length(t)
    if (i==length(t))
        ind_c = find(x==t(i)):length(x);
    else

        ind_c = find(x==t(i)):(find(x==t(i+1))-1);
    end
    if (r(i) == 1)
        ind_pos = [ind_pos,ind_c];
    else if (r(i) == -1)
            ind_neg = [ind_neg,ind_c];
    else
        ind_out = [ind_out,ind_c];
    end
    end

end

%Glucose measurement for increasing value of time
figure
plot(x,y,'Color',[0.7608    0.7059    0.4941],'LineWidth',2);
hold on
sz = 6;
scatter(x(ind_pos),y(ind_pos),sz,'MarkerFaceColor',[0.7098    0.3765    0.3765],'MarkerEdgeColor',[0.7098    0.3765    0.3765]);
scatter(x(ind_neg),y(ind_neg),sz,'MarkerFaceColor',[0.2000    0.4118    0.3686],'MarkerEdgeColor',[0.2000    0.4118    0.3686]);


legend('glucose measurement','increasing trend','decreasing trend','FontSize',10)

end