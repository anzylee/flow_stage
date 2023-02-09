import os
import pandas as pd
from mannings_Qfromh_upstream import *

site = 'sfe_24'
xsect_type = 'bf_to_bf'
path_up_xsect = './'+site+'/'+xsect_type+'/xsect_flow_stage.shp'
path_terrain = './'+site+'/grid/'+site+'.asc'
path_xlsx = './'+site+'/'+xsect_type+'/Q_Stage.xlsx'
initial_water_depth = 1

params = pd.read_excel('./parameters_sfe.xlsx')
for ii in range(0, params.site.__len__() + 1):  # range(0,params.site.__len__()+1)
    if params.site[ii] == site:
        ind = ii
        break
n = params.n[ind]
S0 = params.Slope[ind]

Execute_StackProfile_3d = 1
figure_xsect = 1
unit = 'SI'

Q_array, h_array, A, P, R = mannings_Qfromh_upstream(site, xsect_type, path_up_xsect, path_terrain, initial_water_depth,
                                      n, S0, Execute_StackProfile_3d, figure_xsect, unit)

tmp = np.transpose(np.array([Q_array, h_array]))
df = pd.DataFrame(data=tmp, columns=['Q','h'])
df.to_excel(path_xlsx)