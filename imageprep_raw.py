from pyraf import iraf
import subprocess

def imageprep_raw(zooid, band, i_first, shifts_file, fill):
    coords_file = zooid+'_'+band+'coords.txt'
    flat_file = 'Flat_'+band+'c_sub_test.fits'


    initial_filename_list = []
    initial_filename_list_string = ''

    for index in range(i_first, i_first+5):
        initial_filename_list.append('irr_'+(str(index).zfill(fill))+'_001.fits')
        initial_filename_list_string = initial_filename_list_string + initial_filename_list[-1] + ','

    initial_filename_list_string = initial_filename_list_string[:-1]

    iraf.images()

    skyfile = zooid+'_sky'+band+'_raw.fits'

    iraf.imcombine(initial_filename_list_string, skyfile, combine='median', reject='avsigclip')

    skysub_filename_list = []
    skysub_filename_list_string = ''

    for name in initial_filename_list:
        skysub_filename = name[:-5]+'_skysub_raw.fits'
        iraf.imarith(name,'-',skyfile,skysub_filename)
        skysub_filename_list.append(skysub_filename)
        skysub_filename_list_string = skysub_filename_list_string+skysub_filename+','

    skysub_filename_list_string = skysub_filename_list_string[:-1]

    flattened_filename_list = []
    flattened_filename_list_string = ''

    for name in skysub_filename_list:
        flattened_filename_list.append(name[:-8]+'flattened.fits')
        flattened_filename_list_string = flattened_filename_list_string + flattened_filename_list[-1]+','
        command = 'cp ./'+name+' ./'+flattened_filename_list[-1]
        subprocess.call(command, shell = True)

    flattened_filename_list_string = flattened_filename_list_string[:-1]

    iraf.noao.imred()
    iraf.noao.imred.ccdred()
    iraf.noao.imred.ccdred.ccdproc(flattened_filename_list_string, ccdtype='', fixpix='no', overscan='no', trim ='no', zerocor='no',darkcor='no',flatcor='yes',illumcor='no',fringecor='no',readcor='no',scancor='no',flat=flat_file)

    output_align_filename_list = []
    output_align_filename_list_string = ''

    for i in range(5);
        output_align_filename_list.append(zooid+'_'+band+str(i+1)+'.fits')
        output_align_filename_list_string = output_align_filename_list_string+output_align_filename_list[-1]+','

    output_align_filename_list_string = output_align_filename_list_string[:-1]

    iraf.imalign(flattened_filename_list_string, skysub_filename_list[0], coords_file, output_align_filename_list_string, shifts=shifts_file)

    iraf.imcombine(output_align_filename_list_string, zooid+'_'+band+'_skysubtracted_median.fits',combine='median',reject='none')
    iraf.imcombine(output_align_filename_list_string, zooid+'_'+band+'_skysubtracted_sum.fits',combine='sum',reject='none')
    iraf.imcombine(output_align_filename_list_string, zooid+'_'+band+'_skysubtracted_average.fits',combine='average',reject='none')

    return
