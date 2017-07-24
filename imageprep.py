from pyraf import iraf

def imageprep(zooid, band, shift_file, i_first, indexfill):
    coords_file = zooid+'_'+band+'coords.txt'

    initial_filename_list = []
    initial_filename_list_string = ''

    for index in range(i_first, i_first+5):
        initial_filename_list.append('irr_'+(str(index).zfill(indexfill))+'_001.fits')
        initial_filename_list_string = initial_filename_list_string + initial_filename_list[-1] + ','

    initial_filename_list_string = initial_filename_list_string[:-1]

    iraf.images()

    skyfile = zooid+'_sky'+band+'.fits'

    iraf.imcombine(initial_filename_list_string,skyfile,combine='median',reject='avsigclip')

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
    
    iraf.imcombine(output_align_filename_list_string, zooid+'_'+band+'_skysubtracted_median.fits',combine='median',reject='none')
    iraf.imcombine(output_align_filename_list_string, zooid+'_'+band+'_skysubtracted_sum.fits',combine='sum',reject='none')
    iraf.imcombine(output_align_filename_list_string, zooid+'_'+band+'_skysubtracted_average.fits',combine='average',reject='none')

    return
    
    

