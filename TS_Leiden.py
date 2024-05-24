#This demo runs TS-Leiden algorithm. For details, see:
#
# A. Tastan, C. Escorihuela-Altaba, J. Garcia-Tirado and K. Riesen,
# "Clustering Time Series Data for Personalized Type 1 Diabetes
# Management", in Proc. Int. Conf. Pattern Recognit. Artif. Intell.
# (ICPRAI2024), 2024.
#
# Copyright (C) 2023 Aylin Tastan. All rights reserved.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#
#
# version : 24.05.2024
# author : Aylin Tastan
#######################################################################
import pandas as pd
import numpy as np
from datetime import time
from scipy import stats
#from visibility_graph import visibility_graph
from ts2vg import NaturalVG
from scipy.linalg import issymmetric
import time
import warnings

warnings.filterwarnings("ignore")

def tic():
    #Python version of matlab tic and toc functions
    global startTime_for_tictoc
    startTime_for_tictoc = time.time()

def toc():
    import time
    print (time.time() - startTime_for_tictoc)

def call_data(dir):
    df_cgm = pd.read_csv(f'{dir}/cgm.csv')
    df_basal = pd.read_csv(f'{dir}/basal.csv')
    df_bolus = pd.read_csv(f'{dir}/bolus.csv')
    df_meals = pd.read_csv(f'{dir}/meal.csv')

    return [df_cgm,df_basal,df_bolus,df_meals]

def time_series(data):
    y = data[' Reading']
    datetime_info = pd.to_datetime(data[' Reading taken at'])
    #data['time'] = x.dt.time
    diff_vec = np.zeros(len(datetime_info))
    for i in range (len(datetime_info)-1):
        diff_vec[i+1] = ((datetime_info[i+1] - datetime_info[i]).total_seconds())/60

    x = (datetime_info[0]).minute + np.cumsum(diff_vec)
    return [x,y]

def wTuk(absx,c):
    temp = np.power((absx/c),np.full((absx/c).shape, 2))
    vec_mul = np.ones(absx.shape)
    vec_mul[absx>c] = 0
    wx = np.multiply(np.power(1-temp,np.full((1-temp).shape, 2)),vec_mul)

    return wx

def fast_and_local_moving_of_vertices(A,x,N,p,r,q,c_Tuk,T_weights):
    #Compute modularity and node to cluster edge weights
    TOL = 1e-10
    d = A.sum()
    B = (A-((A.sum(axis=1))*np.transpose([(A.sum(axis=0))]))/d)/d
    if not (issymmetric(B)):
        B = (B+np.transpose(B))/2

    num_clusters = max(p)
    H = np.zeros((num_clusters, num_clusters))

    for i in range(num_clusters):
        idx_clust = np.isin(p,i+1)
        H[:,i] = (B[:,idx_clust]).sum(axis=1)

    while q:

        #Determine the vertex to visit
        v_visit = q[0]
        q.pop(0)
        c_visit = p[v_visit]
        r_visit = r[v_visit]

        #Determine neighbors of the visited vertex robustly
        neigh_v_visit = np.nonzero(A[v_visit,:])
        neigh_v_visit = neigh_v_visit[0]

        if v_visit in neigh_v_visit:
            neigh_v_visit = neigh_v_visit[~np.isin(neigh_v_visit,v_visit)]

        weights = wTuk(abs(x[neigh_v_visit]-np.median(x[neigh_v_visit]))/(1.4826*stats.median_abs_deviation(x[neigh_v_visit])),c_Tuk) #Perform Tukey's weighting function
        neigh_v_visit = neigh_v_visit[weights>T_weights]

        #Calculate the change in quality for neighbors of the visited node
        change_in_quality = H[v_visit,neigh_v_visit]-H[v_visit,c_visit-1]+B[v_visit,v_visit]

        #Determine the local move which provides the maximum quality improvement and perform local move if it provides positive quality improvement
        if (len(change_in_quality)!=0):
            max_quality_change = change_in_quality.max()

            if (max_quality_change>0):
                p_move = p[neigh_v_visit[change_in_quality==max_quality_change]]
                p[v_visit] = p_move[0]
                r_move = r[neigh_v_visit[change_in_quality==max_quality_change]]
                r[v_visit] = r_move[0]

    return[p,r]


def unique_stable(v):

    v_unique = []
    v_refined = []
    ind_clusters = []
    i = 0
    iter = 0
    for item in v:

        if item not in v_unique:
            v_unique.append(item)
            i = i + 1
            ind_clusters.append(iter)
            v_refined.append(i)
        else:
            v_refined.append(v_unique.index(item)+1)

        iter = iter + 1

    return[np.array(v_unique),np.array(ind_clusters),np.array(v_refined)]



def split_vertices(y,p_refined,r_refined,num_clusters):

    ind_outliers = []

    #Split clusters to obtain sequential points associated with the same trend
    for i in range(max(p_refined)):

        ind_cluster_i = np.where(p_refined==i+1)
        ind_cluster_i = ind_cluster_i[0]

        if (len(ind_cluster_i)>1): #if not a singleton cluster
            #print(y[ind_cluster_i])
            r_clust = np.sign(np.diff(y[ind_cluster_i]))

            #print(r_clust)
            r_clust = np.append(r_clust,r_clust[-1])

            #Refine vertices associated with increasing trends
            inc_trend = ind_cluster_i[r_clust>=0]

            if len(inc_trend):
                inc_len = len(inc_trend)
                seq_inc = [np.append(np.append(inc_len,np.diff(inc_trend)),inc_len)==1]
                seq_inc = seq_inc[0]*1
                strseq_inc = ''.join(str(s) for s in seq_inc)
                inc_seq = inc_trend[strseq_inc.find('01'):(strseq_inc.find('10')+1)]

                if len(inc_seq)!=0:
                    r_refined[inc_seq] = 1
            else:
                inc_seq = []


            #Refine vertices associated with decreasing trends
            dec_trend = ind_cluster_i[r_clust<0]

            if len(dec_trend):
                dec_len = len(dec_trend)
                seq_dec = [np.append(np.append(dec_len,np.diff(dec_trend)),dec_len)==1]
                seq_dec = seq_dec[0] * 1
                strseq_dec = ''.join(str(s) for s in seq_dec)
                dec_seq = dec_trend[strseq_dec.find('01'):(strseq_dec.find('10') + 1)]

                if len(dec_seq)!=0:
                    num_clusters = num_clusters + 1
                    p_refined[dec_seq] = num_clusters
                    r_refined[dec_seq] = -1

            else:
                dec_seq = []


            ind_outliers = np.append(ind_outliers,np.setdiff1d(ind_cluster_i,np.append(inc_seq,dec_seq)))

    return[p_refined,r_refined,ind_outliers,num_clusters]


def refine_outliers(y,p_refined,r_refined,ind_outliers,num_clusters):

    trend_act = np.sign(np.diff(np.append(y,y[-1])))
    ind_outliers = np.sort(ind_outliers.astype(int))
    len_p_refined = len(p_refined)
    num_outliers = len(ind_outliers)


    if num_outliers!=0:

        #Define trends associated with the determined outliers
        r_refined[ind_outliers] = trend_act[ind_outliers]

        #Assign first outlier to a cluster
        if (ind_outliers[0]==0):

            if(abs(r_refined[0]-r_refined[1])<=1):
                p_refined[0] = p_refined[1]
            else:
                num_clusters = num_clusters + 1
                p_refined[0] = num_clusters

        #Assign last outlier to a cluster
        if(ind_outliers[-1] == (len_p_refined-1)):

            if(abs(r_refined[len_p_refined-1] - r_refined[len_p_refined-2])<=1):
                p_refined[len_p_refined-1] = p_refined[len_p_refined-2]
            else:
                num_clusters = num_clusters + 1
                p_refined[len_p_refined-1] = num_clusters

        #Assign remaining outliers
        for i in range (1,len(ind_outliers)-1):

            if(abs(r_refined[ind_outliers[i]] - r_refined[ind_outliers[i]-1]) <= 1):
                p_refined[ind_outliers[i]] = p_refined[ind_outliers[i]-1]

            elif(abs(r_refined[ind_outliers[i]] - r_refined[ind_outliers[i]+1])<=1):
                p_refined[ind_outliers[i]] = p_refined[ind_outliers[i]+1]

            else:
                num_clusters = num_clusters + 1
                p_refined[ind_outliers[i]] = num_clusters


    return[p_refined,r_refined,num_outliers]


def refine_partition(y,p,r):

    #Determine the number of clusters and order cluster assignments
    [p_unique,_,p_refined] = unique_stable(p)
    num_clusters = len(p_unique)

    #Split vertices based on their trends
    [p_refined, r_refined, ind_outliers, num_clusters] = split_vertices(y,p_refined,r,num_clusters)

    if (len(ind_outliers)>0):
        #Refine information corresponding to outliers
        [p_refined,r_refined,num_outliers] = refine_outliers(y.array,p_refined,r_refined,ind_outliers,num_clusters)
    else:
        num_outliers = 0


    return[p_refined,r_refined,num_outliers]


def aggregate_graph(A,y,x,p_refined,r_refined,max_neigh_num):

    #Obtain the number of clusters and their indices
    [_,ind_clusters, p_refined] = unique_stable(p_refined)
    N_agg = max(p_refined)

    #Aggregate all variables
    x_agg = x[ind_clusters]
    r_agg = r_refined[ind_clusters]
    y_agg = y[ind_clusters]
    y_agg.reset_index(inplace=True, drop=True)
    p_agg = p_refined[ind_clusters]

    #Aggregate the affinity matrix
    ind_clusters = np.append(ind_clusters,len(p_refined))
    A_agg = np.zeros((N_agg,N_agg))


    for i in range(N_agg):
       for j in range(i,N_agg):

           if(j > i+max_neigh_num):
               break
           else:
               A_agg[i,j] = (A[ind_clusters[i]:(ind_clusters[j+1]),ind_clusters[j]:(ind_clusters[j+1])]).sum()


    A_agg = A_agg + np.tril(np.transpose(A_agg),-1)


    return[A_agg,y_agg,x_agg,N_agg,p_agg,r_agg]


#Perform TS Leiden for risky time determination
#Define directory and obtain data
dir = '/home/atastan/Desktop/RESEARCH_PRG/UDEM_Collaboration/DataCollection_DCLP6/02f March_1_2022_Data_Collection/85201_20220118_20220216'
[df_cgm,df_basal,df_bolus,df_meals] = call_data(dir)
#print(df_cgm.columns)

#Determine the time series which will be analyzed
[x,y] = time_series(df_cgm)
x_ini = x
y_ini = y

#Step 0: Initialization
g = NaturalVG(directed=None).build(y)
A = g.adjacency_matrix() #adjacency matrix for natural visibility graph
N = len(y)
p = np.array(range(1,N+1))  #initial partition vector
r = np.zeros(N)  #initial trend vector
q = list(range(N))  #vertex indices to visit
c_Tuk = 3.4447 #Tukey's tuning constant for real valued data
T_weights = 0 #weight threshold (default is zero and it means that the threshold is not used)
max_neigh_num = 10 #user-defined nearest neighbor value for graph aggregation
TOL = 5 #tolerence level for the number of outliers (alternatively number of iterations can be determined)


while 1:

    #Step 1: Fast and Robust Local Moving of Vertices
    [p,r] = fast_and_local_moving_of_vertices(A,x,N,p,r,q,c_Tuk,T_weights)

    #Step 2: Refinement of the Partition
    [p,r,num_outliers] = refine_partition(y,p,r)

    #Step 3: Aggregation of the Graph based on the Refined Partition
    [A,y,x,N,p,r] = aggregate_graph(A,y,x,p,r,max_neigh_num)

    q = np.where(r==0)
    q = (q[0]).tolist()

    if not (len(q)>TOL):
        break


#Summary of clustering and trend determination results
p_res = np.zeros(len(x_ini))
r_res = np.zeros(len(x_ini))
x_res = np.append(x,x_ini[-1]+1)
for i in range(len(p)):
    ind_clust_i = np.where(np.logical_and(x_ini>=x_res[i], x_ini<x_res[i+1]))
    p_res[ind_clust_i] = i+1
    r_res[ind_clust_i] = r[i]

res_mat = (np.stack([p_res,r_res], axis=1)).astype(int)










