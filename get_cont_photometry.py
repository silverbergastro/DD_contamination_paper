import numpy as np
from astropy.io import ascii
from astropy.io import table

def get_cont_photometry(zooid, band, num, targ_x, targ_y):
    #Get name of file from zooid, band, which version to use
    filename = zooid+'_'+band+'final.fits.mag.'+str(num)
    data = ascii.read(filename)    #Read in data

    target_x_est = targ_x
    target_y_est = targ_y

    length = data['ID'][-1]    #Get number of objects in file

    for i in range(length):
        #Identify the target in the file
        if (abs(data['YCENTER'][i] - target_y_est) < 0.1) & (abs(data['XCENTER'][i] - target_x_est) < 0.1):
            target_id = data['ID'][i]
            target_x = data['XCENTER'][i]
            target_y = data['YCENTER'][i]
            target_radius = np.sqrt(data['AREA'][i]/np.pi)

    #Create table for contaminants
    cont_table = data[:0].copy()
    cont_too_close_table = data[:0].copy()

    #Create table for contaminants with overlap
    cont_table.add_row(data[target_id-1])
    cont_too_close_table.add_row(data[target_id-1])

    #Add stars to appropriate tables, skipping the target
    for i in range(length):
        if data['ID'][i] == target_id:
            continue
        else:
            xcen = data['XCENTER'][i]
            ycen = data['YCENTER'][i]
            xdist = xcen - target_x
            ydist = ycen - target_y
            dist = np.sqrt(xdist**2 + ydist**2)
            if dist < 60.:
                cont_table.add_row(data_test[i])
            cont_radius = np.sqrt(data['AREA'][i]/np.pi)
            if dist < (cont_radius + target_radius):
                cont_too_close_table.addrow(data[i])

    cont_filename = zooid+'_'+band+'final_contaminants.dat'
    cont_prox_filename = zooid+'_'+band+'final_contaminants_too_close.dat'

    #Write to file
    if np.array(cont_table['IMAGE']).size > 0:
        ascii.write(cont_too_close_table, cont_prox_filename)
    else:
        print "No contaminants"

    if np.array(cont_too_close_table['IMAGE']).size >0:
        ascii.write(cont_table, cont_filename)
    else:
        print "No overlapping contaminants"

    return
