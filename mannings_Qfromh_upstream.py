import os
import simpledbf
import arcpy
import _thread
import matplotlib.pyplot as plt
import numpy as np


# This script is from F:\tuflow-SC\py_modules

def mannings_Qfromh_upstream(site, xsect_type, path_up_xsect, path_terrain, water_depth,
                             n, S0, Execute_StackProfile_3d, figure_xsect, unit):
    # This python module calculates the flow discharge using the water depth at a certain cross-section
    # The required inputs are:
    #   path_up_xsect         - path to the transect line (XX.shp, shape file)
    #   path_terrain            - path to the terraian (XX.asc)
    #   water_stage             - the water stage at the xsect interested (e.g. 1st riffle-crest)
    #   n                       - manning's n
    #   S0                      - Slope of the channel
    #   Execute_StackProfile_3d - 1 if you need to generate xsect profile in xlsx
    #   figure_xsect            - 1 if you want to see the xsect profile and water stage
    #
    # The output is:
    #   Q                       - Calculated discharge for the site
    #   A                       - Flow area when the discharge = Q
    #   P                       - Perimeter when the discharge = Q
    #   R                       - Hydraulic radius when the discharge = Q

    Q_array, h_array = [], []
    # path to upstream xsect shp file
    xsectshp1 = path_up_xsect

    print(xsectshp1)

    # Load Raster DEM
    terrain = path_terrain

    # Define projection
    dsc = arcpy.Describe(xsectshp1)
    coord_sys = dsc.spatialReference
    #arcpy.DefineProjection_management(terrain, coord_sys)

    # Stack Profile
    # upstream: at the first riffle-crest
    xsecttab1 = './'+site+'/'+xsect_type+'/xsect_table.dbf'
    print(xsecttab1)

    if Execute_StackProfile_3d:
        if not (os.path.isfile(xsectshp1)):
            print(xsectshp1 + ' is not found. Please create a polyline in arcGIS and save it in this location.')
            _thread.interrupt_main()
        #    _thread.interrupt_main()

        if os.path.isfile(xsecttab1):
            os.remove(xsecttab1)

        # Execute Stack Profile
        arcpy.CheckOutExtension("3D")
        arcpy.StackProfile_3d(xsectshp1, profile_targets=[terrain], out_table=xsecttab1)

    xsectdbf1 = simpledbf.Dbf5(xsecttab1)
    xsectdfst1 = xsectdbf1.to_dataframe()
    xsectdf1 = xsectdfst1
    #Line_IDs = range(0, max(xsectdf1['LINE_ID'])+1)
    Line_IDs = xsectdf1['LINE_ID'].unique()

    for Line_ID in Line_IDs:
        print("Line ID = ", str(Line_ID))
        # Construct a functional relationship between A and h
        x = np.array(xsectdf1.loc[xsectdf1['LINE_ID']==Line_ID]['FIRST_DIST'])
        z = np.array(xsectdf1.loc[xsectdf1['LINE_ID']==Line_ID]['FIRST_Z'])
        elevation = min(z)

        water_stage = elevation + water_depth

        z0 = z - water_stage
        #print(z0)
        ind = []

        if figure_xsect == 1:

            # Figure, at the first riffle-crest
            plt.figure(1)
            plt.plot(x, z, '-')
            plt.plot([np.min(x), np.max(x)], [water_stage, water_stage], '-')
            plt.grid()
            plt.xlabel('Lateral Distance '+'('+unit+')')
            plt.ylabel('Elevation '+'('+unit+')')
            plt.title('Cross-sectional profile at the upstream (or downstream)')
            #plt.show()
            ini_path_fig = './' + site + '/' + xsect_type + '/xsect_profile_' + str(Line_ID)
            plt.savefig(ini_path_fig)
            plt.close()

        new_water_stage = float(input("Enter your water stage (check "+ini_path_fig+"): "))

        water_stage = new_water_stage
        z0 = z - water_stage

        # Calculate Flow discharge using Manning's n

        for ii in range(0,z.__len__()-1):
            if np.sign(z0[ii]*z0[ii+1]) < 0:
                ind.append(ii)
        A = 0
        P = 0

        for ii in range(0, ind.__len__(), 2):

            m1 = (z0[ind[ii]] - z0[ind[ii] + 1]) / (x[ind[ii]] - x[ind[ii] + 1])
            xi1 = (-z0[ind[ii]] + m1 * x[ind[ii]]) / m1

            m2 = (z0[ind[ii+1]] - z0[ind[ii+1] + 1]) / (x[ind[ii+1]] - x[ind[ii+1] + 1])
            xi2 = (-z0[ind[ii+1]] + m2 * x[ind[ii+1]]) / m2

            X = np.hstack((xi1, x[ind[ii] :ind[ii+1]+1], xi2))
            Z = -np.hstack((0, z0[ind[ii] :ind[ii+1]+1], 0))
            dA = np.trapz(Z, x=X)

            dx = X[1:] - X[:-1]
            dz = Z[1:] - Z[:-1]
            dP = np.sum((dx ** 2 + dz ** 2) ** (1 / 2))
            print('dA='+str(dA))
            A = A + dA
            P = P + dP


        R = A/P
        if unit == 'feet':
            Q = (1.49/n)*A*R**(2/3)*(S0)**(1/2)
        elif unit in ['SI', 'meter']:
            Q = (1 / n) * A * R ** (2 / 3) * (S0) ** (1 / 2)

        if figure_xsect == 1:

            # Figure, at the first riffle-crest
            plt.figure(2)
            plt.plot(x, z, '-')
            plt.plot([np.min(x), np.max(x)], [water_stage, water_stage], '-')
            plt.grid()
            plt.xlabel('Lateral Distance '+'('+unit+')')
            plt.ylabel('Elevation '+'('+unit+')')
            plt.title('Cross-sectional profile at the upstream (or downstream)')
            #plt.show()
            plt.savefig('./' + site + '/' + xsect_type + '/xsect_profile_final_' + str(Line_ID))
            plt.close()

        print('#### Q = '+str(Q)+' ####')
        print('#### h = ' + str(water_stage) + ' ####')
        Q_array = np.append(Q_array, Q)
        h_array = np.append(h_array, water_stage)

    return Q_array, h_array, A, P, R
