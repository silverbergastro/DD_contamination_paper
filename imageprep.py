from pyraf import iraf
import numpy as np
import time
from astropy.io import fits

def imageprep(zooid, band, shift_file, i_first, indexfill):
    coords_file = zooid+'_'+band+'coords.txt'
    iraf.images()

    print zooid, band, time.time()

    #initial_filename_list_sky = []
    #initial_filename_list_sky_string = ''

    initial_filename_list = []
    initial_filename_list_string = ''

    for index in range(i_first, i_first+5):
        initial_filename_list.append('irr_'+(str(index).zfill(indexfill))+'_001.fits')
        initial_filename_list_string = initial_filename_list_string + initial_filename_list[-1] + ','
        #if index >= i_first_science:
        #    initial_filename_list.append(initial_filename_list_sky[-1])
        #    initial_filename_list_string = initial_filename_list_string+initial_filename_list[-1]+','

    #initial_filename_list_sky_string = initial_filename_list_sky_string[:-1]
    initial_filename_list_string = initial_filename_list_string[:-1]

    file_array_list = []
    for filename in initial_filename_list:
        hdulist = fits.open(filename)
        file_array_list.append(hdulist[0].data)
        hdulist.close()

    all_file_median = np.median(file_array_list)

    individual_file_median = [np.median(i) for i in file_array_list]
    individual_file_median_of_medians = np.median(individual_file_median)

    scaling_used = False

    if (np.max(individual_file_median) - np.min(individual_file_median)) > 0.1*all_file_median:
        scaling_used = True
        print "Large variance in images--scaling required"
        scaled_filename_list = []
        scaled_filename_list_string = ''
        for i in range(5):
            filename_init = initial_filename_list[i]
            filename_init_median = individual_file_median[i]

            scaled_filename_list.append(filename_init[:-5]+'_scaled.fits')
            scaled_filename_list_string = scaled_filename_list_string + scaled_filename_list[-1]+','

            outlier_test = np.abs(filename_init_median - all_file_median)
            new_margin = filename_init_median - all_file_median
            if outlier_test > 100.:
                while np.abs(new_margin) > 100.:
                    new_margin = 0.1*new_margin
                new_median = all_file_median + new_margin
                scale_factor = new_median/filename_init_median
                pedestal_add = new_median - filename_init_median
            else:
                scale_factor = 1.
                pedestal_add = 1.

            iraf.imarith(filename_init,'*',scale_factor,scaled_filename_list[-1])

        scaled_filename_list_string = scaled_filename_list_string[:-1]
            

    #skyfile = zooid+'_sky'+band+'_'+str(hsigma)+'_'+str(lsigma)+'_'+reject_type+'.fits'
    skyfile = zooid+'_sky'+band+'.fits'

    #input_hsigma = str(hsigma)
    #input_lsigma = str(lsigma)

    #iraf.imcombine(initial_filename_list_string,skyfile,combine='median',reject=reject_type,lsigma=input_lsigma,hsigma=input_hsigma)
    if scaling_used:
        iraf.imcombine(scaled_filename_list_string,skyfile,combine='median',reject='avsigclip',hsigma='3.',rdnoise='2.',gain='25.')
    else:
        iraf.imcombine(initial_filename_list_string, skyfile, combine='median', reject='avsigclip',hsigma='3.',rdnoise='2.',gain='25.')

    skysub_filename_list = []
    skysub_filename_list_string = ''

    for name in initial_filename_list:
        skysub_filename = name[:-5]+'_skysub.fits'
        print skysub_filename
        iraf.imarith(name,'-',skyfile,skysub_filename)
        skysub_filename_list.append(skysub_filename)
        skysub_filename_list_string = skysub_filename_list_string+skysub_filename_list[-1]+','

    skysub_filename_list_string = skysub_filename_list_string[:-1]

    output_align_filename_list = []
    output_align_filename_list_string = ''

    for i in range(5):
        output_align_filename_list.append(zooid+'_'+band+str(i+1)+'.fits')
        output_align_filename_list_string = output_align_filename_list_string + output_align_filename_list[-1] + ','

    output_align_filename_list_string = output_align_filename_list_string[:-1]

    print skysub_filename_list[0]

    iraf.imalign(skysub_filename_list_string, skysub_filename_list[0], coords_file, output_align_filename_list_string,shifts=shift_file)
    
    iraf.imcombine(output_align_filename_list_string, zooid+'_'+band+'_skysubtracted_median.fits',combine='median',reject='crreject')
    iraf.imcombine(output_align_filename_list_string, zooid+'_'+band+'_skysubtracted_sum.fits',combine='sum',reject='crreject')
    iraf.imcombine(output_align_filename_list_string, zooid+'_'+band+'_skysubtracted_average.fits',combine='average',reject='crreject')

    return
