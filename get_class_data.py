import numpy as np
import json
import matplotlib.pyplot as plt
from pprint import pprint
import re
import matplotlib as mpl
import time

class Subject:
    def __init__(self, zooniverse_id, wise_id, subject_dict, classlist):
        self.zooniverse_id = zooniverse_id
        self.wise_id = wise_id
        self.num_classifiers = len(classlist)
        self.classification_dict = {}
        self.classification_dict[u'multi'] = 0
        self.classification_dict[u'good'] = 0
        self.classification_dict[u'oval'] = 0
        self.classification_dict[u'empty'] = 0
        self.classification_dict[u'extended'] = 0
        self.classification_dict[u'shift'] = 0
        self.state = subject_dict[u'state']
        self.jmag = None
        self.jmagerr = None
        self.hmag = None
        self.hmagerr = None
        self.kmag = None
        self.kmagerr = None
        self.w1mag = None
        self.w1magerr = None
        self.w2mag = None
        self.w2magerr = None
        self.w3mag = None
        self.w3magerr = None
        self.w4mag = None
        self.w4magerr = None

        self.eq_coords = None
        self.gal_coords = None

        #self.wise_id = subject_dict[u'metadata'][u'wise_id']
        if len(subject_dict[u'coords']) > 0:
            self.eq_coords = subject_dict[u'coords']
        else:
            get_coords_string = subject_dict[u'metadata'][u'wise_id']
            self.eq_coords = get_eq_coords(get_coords_string)
        #self.eq_coords = subject_dict[u'coords']
        self.gal_coords = get_gal_coords(self.eq_coords)
        #self.jmag = jmag
        #self.jmagerr = jmagerr

        #self.eq_coords = None
        #self.gal_coords = None
        self.majmult = False
        self.majgood = False
        self.other = False
        #if ((self.classification_dict[u'good']/self.num_classifiers) > 0.5):
        #    self.majgood = True
        #elif ((self.classification_dict[u'multi']/self.num_classifiers) > 0.5):
        #    self.majmult = True
        #else:
        #    self.other = True
        self.subject_dict = subject_dict
        self.classlist = classlist

    def __str__(self):
        printlist = []
        printlist.append(str(self.zooniverse_id))
	printlist.append(str(self.wise_id))
        #printlist.append('['+str(self.eq_coords[0])+','+str(self.eq_coords[1])+']')
        #printlist.append('['+str(self.gal_coords[0])+','+str(self.gal_coords[1])+']') 
        printlist.append(str(self.num_classifiers))
        printlist.append(str(self.classification_dict[u'good']))
        printlist.append(str(self.classification_dict[u'multi']))
        printlist.append(str(self.classification_dict[u'oval']))
        printlist.append(str(self.classification_dict[u'empty']))
        printlist.append(str(self.classification_dict[u'extended']))
        printlist.append(str(self.classification_dict[u'shift']))
        printlist.append(str(self.state))
        #printlist.append(str(self.jmag))
        #printlist.append(str(self.jmagerr))
        #printlist.append(str(self.hmag))
        #printlist.append(str(self.hmagerr))
        #printlist.append(str(self.kmag))
        #printlist.append(str(self.kmagerr))
        #printlist.append(str(self.w1mag))
        #printlist.append(str(self.w1magerr))
        #printlist.append(str(self.w2mag))
        #printlist.append(str(self.w2magerr))
        #printlist.append(str(self.w3mag))
        #printlist.append(str(self.w3magerr))
        #printlist.append(str(self.w4mag))
        #printlist.append(str(self.w4magerr))
   
        #printlist.append("jmag : " + str(self.jmag) + " +/m " + str(self.jmagerr))
        #if self.majmult is True:
        #    printlist.append("Maj Mult")
        #elif self.majgood is True:
        #    printlist.append("Maj Good")
        #else:
        #    printlist.append("Other")

        s = ""
        for entry in printlist:
            s = s + entry + ","

        s=s[:-1]

        return s

def driver(filename1, filename2,  directory):
    start_time = time.time()
    print start_time
    print "Importing class data"
    class_data = import_class_data(filename1)
    subject_data = import_data(filename2)
    #phot_data = import_phot_data(filename3)

    #photometry_mags_data, photometry_errs_data, eq_coords_dict, gal_coords_dict = import_photometry_data(filename3)

    print "classifications :", len(class_data)
    print "subjects :", len(subject_data)

    #print "photometry :", len(photometry_mags_data.keys())

    #pprint(subject_data[43042])
    #pprint(subject_data[43043])
    #pprint(subject_data[43044])

    #objects_vec, objects_map = build_objects(subject_data, class_data, photometry_mags_data, photometry_errs_data, eq_coords_dict, gal_coords_dict)
    objects_vec, objects_map_zoo, objects_map_wise = build_objects(subject_data, class_data)
   
    #completed_objects = get_completed_zooids(objects_vec)

    del subject_data
    del class_data

    #goodlist, multilist, otherlist = get_totals(completed_objects)

    #rates, rate_errs = coord_number_histogram(goodlist, multilist, otherlist, directory)

    #doneness = coord_density_histogram(rates, rate_errs, directory)
    #print doneness

    print_to_file(objects_vec, directory)

    #doneness1 = lat_long_heatmap(goodlist, multilist, otherlist, directory)

    print "Done"
    end_time = time.time()
    print end_time
    runtime = end_time - start_time
    print runtime

    print "Done"

    return

#def import_phot_data(filename):
    #datavec = []
#    phot_dict = {}
#    with open(filename) as f:
#        for line in f:
#            load = json.loads(line)
#            #datavec.append(json.loads(line))
#            wise_id = load[u'wiseid']
#            phot_dict[wise_id] = load
#            if len(phot_dict.keys()) < 2:
#                print wise_id
#            if (len(phot_dict.keys()) % 100) < 1:
#                print len(phot_dict.keys())
#    return phot_dict

def import_class_data(filename):
    datavec = []
    bad_ids = [u'AWI0000x7o', u'AWI000035p', u'AWI00006vl', u'AWI0000hlg']
    with open(filename) as f:
        for line in f:
	    load = json.loads(line)
            if not load[u'tutorial']:
                if len(load[u'subjects']) > 0:
                    if load[u'subjects'][0][u'zooniverse_id'] not in bad_ids:
                        datavec.append(json.loads(line))
            if (len(datavec) % 100) <1:
                print len(datavec)
    return datavec

def import_data(filename):
    datavec = []
    bad_ids = [u'AWI0000x7o', u'AWI000035p', u'AWI00006vl', u'AWI0000hlg']
    with open(filename) as f:
        for line in f:
            load = json.loads(line)
            if load[u'zooniverse_id'] in bad_ids:
                print 'Bad ID skipped'
                continue             
            datavec.append(json.loads(line))
            if (len(datavec) % 100) < 1:
                print len(datavec)
    return datavec

def build_objects(subject_vec, classifications_vec):
    objects_vec = []
    objects_map_zoo = {}
    #objects_map_wise = {}
    i = 0

    build_start_time = time.time()
    wise_id_list = []
    
    for entryitem in subject_vec:
        wise_id = entryitem[u'metadata'][u'wise_id']
        wise_id_list.append(wise_id)

    objects_map_wise = {key: [] for key in set(wise_id_list)}

    for entry in subject_vec:
        name = entry[u'zooniverse_id']
        wise_id = entry[u'metadata'][u'wise_id']
        classlist = []
        if (i%100) < 1:
            print i, name
        if ((i%20000) < 1) and (i>0):
            print 'Average subject time', (time.time() - build_start_time)/i
        if name == u'AWI0000x7o':
            pprint(entry)
            print i
        x = Subject(name, wise_id, entry, classlist)
        objects_vec.append(x)
        objects_map_zoo[name] = (len(objects_vec) - 1)
        #if wise_id not in objects_map_wise.keys():
        #    objects_map_wise[wise_id] = []
        objects_map_wise[wise_id].append(len(objects_vec)-1)
        i += 1

    print 'Subjects read into Objects'

    class_start_time = time.time()

    j = 0
    for item in classifications_vec:
        if len(item[u'subjects']) > 0:
            classname = item[u'subjects'][0][u'zooniverse_id']
            #print classname
            subj_index = objects_map_zoo[classname]
            #print len(objects_vec[subj_index].classlist)
            objects_vec[subj_index].classlist.append(item)
            objects_vec[subj_index].num_classifiers += 1
            #print len(objects_vec[subj_index].classlist)
            #print objects_vec[subj_index].num_classifiers
            j += 1
            if (j % 20000) < 1:
                print (time.time() - class_start_time)/j

            subitem = item[u'annotations']
            for entry in subitem:
                if u'classified_as' in entry.keys():
                    res = entry[u'classified_as']
                    objects_vec[subj_index].classification_dict[res] += 1

    print 'Classifications read into Objects'

    #for item in phot_vec:
    #    wise_id = item[u'wiseid']
    #    if wise_id in objects_map_wise.keys():
    #        subj_indices = objects_map_wise[wise_id]
    #        for subj_index in subj_indices:
    #            objects_vec[subj_index].jmag = item[u'jmag']
    #            objects_vec[subj_index].jmagerr = item[u'jmagerr']
    #            objects_vec[subj_index].hmag = item[u'hmag']
    #            objects_vec[subj_index].hmagerr = item[u'hmagerr']
    #            objects_vec[subj_index].kmag = item[u'kmag']
    #            objects_vec[subj_index].kmagerr = item[u'kmagerr']
    #            objects_vec[subj_index].w1mag = item[u'w1mag']
    #            objects_vec[subj_index].w1magerr = item[u'w1magerr']
    #            objects_vec[subj_index].w2mag = item[u'w2mag']
    #            objects_vec[subj_index].w2magerr = item[u'w2magerr']
    #            objects_vec[subj_index].w3mag = item[u'w3mag']
    #            objects_vec[subj_index].w3magerr = item[u'w3magerr']
    #            objects_vec[subj_index].w4mag = item[u'w4mag']
    #            objects_vec[subj_index].w4magerr = item[u'w4magerr']
    #            objects_vec[subj_index].eq_coords = item[u'eq_coords']
    #            objects_vec[subj_index].gal_coords = item[u'gal_coords']

#    ticker = 0
#    with open(filename3) as f:
#        for line in f:
#            ticker += 1
#            load = json.loads(line)
#            #datavec.append(json.loads(line))
#            wise_id = load[u'wiseid']
#            if wise_id in objects_map_wise.keys():
#                subj_indices = objects_map_wise[wise_id]
#	        for subj_index in subj_indices:
#                    objects_vec[subj_index].jmag = load[u'jmag']
#                    objects_vec[subj_index].jmagerr = load[u'jmagerr']
#                    objects_vec[subj_index].hmag = load[u'hmag']
#                    objects_vec[subj_index].hmagerr = load[u'hmagerr']
#                    objects_vec[subj_index].kmag = load[u'kmag']
#                    objects_vec[subj_index].kmagerr = load[u'kmagerr']
#                    objects_vec[subj_index].w1mag = load[u'w1mag']
#                    objects_vec[subj_index].w1magerr = load[u'w1magerr']
#                    objects_vec[subj_index].w2mag = load[u'w2mag']
#                    objects_vec[subj_index].w2magerr = load[u'w2magerr']
#                    objects_vec[subj_index].w3mag = load[u'w3mag']
#                    objects_vec[subj_index].w3magerr = load[u'w3magerr']
#                    objects_vec[subj_index].w4mag = load[u'w4mag']
#                    objects_vec[subj_index].w4magerr = load[u'w4magerr']
#                    objects_vec[subj_index].eq_coords = load[u'eq_coords']
#                    objects_vec[subj_index].gal_coords = load[u'gal_coords']
#            if (ticker%100) < 1:
#                print ticker
 
    #for wise_id in objects_map_wise.keys():
    #    phot = phot_dict[wise_id]
    #    subj_indices = objects_map_wise[wise_id]
    #    for subj_index in subj_indices:
    #        objects_vec[subj_index].jmag = item[u'jmag']
    #        objects_vec[subj_index].jmagerr = item[u'jmagerr']
    #        objects_vec[subj_index].hmag = item[u'hmag']
    #        objects_vec[subj_index].hmagerr = item[u'hmagerr']
    #        objects_vec[subj_index].kmag = item[u'kmag']
    #        objects_vec[subj_index].kmagerr = item[u'kmagerr']
    #        objects_vec[subj_index].w1mag = item[u'w1mag']
    #        objects_vec[subj_index].w1magerr = item[u'w1magerr']
    #        objects_vec[subj_index].w2mag = item[u'w2mag']
    #        objects_vec[subj_index].w2magerr = item[u'w2magerr']
    #        objects_vec[subj_index].w3mag = item[u'w3mag']
    #        objects_vec[subj_index].w3magerr = item[u'w3magerr']
    #        objects_vec[subj_index].w4mag = item[u'w4mag']
    #        objects_vec[subj_index].w4magerr = item[u'w4magerr']
    #        objects_vec[subj_index].eq_coords = item[u'eq_coords']
    #        objects_vec[subj_index].gal_coords = item[u'gal_coords']

    #print 'Photometry read into Objects'
           

    #for designation in mags_dict.keys():
    #    subj_indices = objects_map_wise[designation]
    #    for number in subj_indices:
    #        objects_vec[number].jmag = mags_dict[designation][0]
    #        objects_vec[number].jmagerr = mags_err_dict[designation][1]
    #        objects_vec[number].hmag = mags_dict[designation][2]
    #        objects_vec[number].hmagerr = mags_err_dict[designation][3]
    #        objects_vec[number].kmag = mags_dict[designation][4]
    #        objects_vec[number].kmagerr = mags_err_dict[designation][5]
    #        objects_vec[number].w1mag = mags_dict[designation][6]
    #        objects_vec[number].w1magerr = mags_err_dict[designation][7]
    #        objects_vec[number].w2mag = mags_dict[designation][8]
    #        objects_vec[number].w2magerr = mags_err_dict[designation][9]
    #        objects_vec[number].w3mag = mags_dict[designation][10]
    #        objects_vec[number].w3magerr = mags_err_dict[designation][11]
    #        objects_vec[number].w4mag = mags_dict[designation][12]
    #        objects_vec[number].w4magerr = mags_err_dict[designation][13]
    #        objects_vec[number].eq_coords = eqcoords_dict[designation]
    #        objects_vec[number].gal_coords = galcoords_dict[designation]

    #'''Combine Zoo IDs with the same WISE id'''

    return objects_vec, objects_map_zoo, objects_map_wise


#def build_objects(subject_vec, classifications_vec, mags_dict, errs_dict, eqcoords_dict, galcoords_dict):
#    objects_vec = []
#    objects_map_zoo = {}
#    objects_map_wise = {}
#    i = 0
#    for entry in subject_vec:
#        name = entry[u'zooniverse_id']
#        wise_id = entry[u'metadata'][u'wise_id']
#        classlist = []
#        print name
#        if name == u'AWI0000x7o':
#            pprint(entry)
#            print i
#        x = Subject(name, entry, classlist)
#        objects_vec.append(x)
#        objects_map_zoo[name] = (len(objects_vec) - 1)
#        if wise_id not in objects_map_wise.keys():
#            objects_map_wise[wise_id] = []
#        objects_map_wise[wise_id].append(len(objects_vec)-1)
#        i += 1

#    print 'Subjects read into Objects'

#    for item in classifications_vec:
#        if len(item[u'subjects']) > 0:
#            classname = item[u'subjects'][0][u'zooniverse_id']
            #print classname
#            subj_index = objects_map[classname]
            #print len(objects_vec[subj_index].classlist)
#            objects_vec[subj_index].classlist.append(item)
#            objects_vec[subj_index].num_classifiers += 1
            #print len(objects_vec[subj_index].classlist)
            #print objects_vec[subj_index].num_classifiers

#            subitem = item[u'annotations']
#            for entry in subitem:
#                if u'classified_as' in entry.keys():
#                    res = entry[u'classified_as']
#                    objects_vec[subj_index].classification_dict[res] += 1

#    print 'Classifications read into Objects'

#    for designation in mags_dict.keys():
#        subj_indices = objects_map_wise[designation]
#        for number in subj_indices:
#            objects_vec[number].jmag = mags_dict[designation][0]
#            objects_vec[number].jmagerr = mags_err_dict[designation][1]
#            objects_vec[number].hmag = mags_dict[designation][2]
#            objects_vec[number].hmagerr = mags_err_dict[designation][3]
#            objects_vec[number].kmag = mags_dict[designation][4]
#            objects_vec[number].kmagerr = mags_err_dict[designation][5]
#            objects_vec[number].w1mag = mags_dict[designation][6]
#            objects_vec[number].w1magerr = mags_err_dict[designation][7]
#            objects_vec[number].w2mag = mags_dict[designation][8]
#            objects_vec[number].w2magerr = mags_err_dict[designation][9]
#            objects_vec[number].w3mag = mags_dict[designation][10]
#            objects_vec[number].w3magerr = mags_err_dict[designation][11]
#            objects_vec[number].w4mag = mags_dict[designation][12]
#            objects_vec[number].w4magerr = mags_err_dict[designation][13]
#            objects_vec[number].eq_coords = eqcoords_dict[designation]
#            objects_vec[number].gal_coords = galcoords_dict[designation]'''

#    '''Combine Zoo IDs with the same WISE id'''

#    '''for designation in mags_dict.keys():
#        if len(objects_map_wise[designation]) > 1:
#           indices = objects_map_wise[designation]
#           for index in indices[1:]:
#               objects_vec[indices[0]].classification_dict[u'multi'] += objects_vec[index].classification_dict[u'multi']
#               objects_vec[indices[0]].classification_dict[u'good'] += objects_vec[index].classification_dict[u'good']
#               objects_vec[indices[0]].classification_dict[u'oval'] += objects_vec[index].classification_dict[u'oval']
#               objects_vec[indices[0]].classification_dict[u'empty'] += objects_vec[index].classification_dict[u'empty']
#               objects_vec[indices[0]].classification_dict[u'extended'] += objects_vec[index].classification_dict[u'extended']
#               objects_vec[indices[0]].classification_dict[u'shift'] += objects_vec[index].classification_dict[u'shift']
#               objects_vec[indices[0]].num_classifiers += objects_vec[index].num_classifiers

#           for index in indices[1:]:
#               removed_zoo_id = objects_vec[index].zooniverse_id
#               del objects_vec[index]
#               del objects_map[removed_zoo_id]
           
#           objects_map_wise[designation] = list(objects_map_wise[designation][0])

#    return objects_vec, objects_map_zoo, objects_map_wise'''

    

def get_completed_zooids(subject_vec):
    completed_zooids = []
    for entry in subject_vec:
        if entry.state == u'complete' and entry.num_classifiers > 0:
            #entry.num_classifiers = len(entry.classlist)
            completed_zooids.append(entry)
    return completed_zooids

def get_gal_coords(eq_coords_vec):
    ra_deg = eq_coords_vec[0]
    dec_deg = eq_coords_vec[1]

    ra_rad = ra_deg*np.pi/180.
    dec_rad = dec_deg*np.pi/180.

    alphag_deg = (12+(51.4/60))*15
    deltag_deg = 27.13

    alphab_deg = (17+(45.6/60))*15
    deltab_deg = -28.94
    
    #convert to radians
    alphag_rad = alphag_deg*np.pi/180
    deltag_rad = deltag_deg*np.pi/180
    alphab_rad = alphab_deg*np.pi/180
    deltab_rad = deltab_deg*np.pi/180

    cosBK = (np.sin(deltab_rad)*np.cos(deltag_rad)) - (np.cos(deltab_rad)*np.sin(deltag_rad)*np.cos(alphag_rad - alphab_rad))
    BK = np.arccos(cosBK)

    sinb = (np.sin(deltag_rad)*np.sin(dec_rad)) + (np.cos(deltag_rad)*np.cos(dec_rad)*np.cos(ra_rad-alphag_rad))
    cosbsinang = np.cos(dec_rad)*np.sin(ra_rad - alphag_rad)
    cosbcosang = (np.cos(deltag_rad)*np.sin(dec_rad)) - (np.sin(deltag_rad)*np.cos(dec_rad)*np.cos(ra_rad-alphag_rad))

    b_rad = np.arcsin(sinb)
    cosb = np.cos(b_rad)

    sinang = cosbsinang/cosb
    cosang = cosbcosang/cosb

    sinang_use = -sinang
    cosang_use = +cosang

    ang1 = np.arcsin(sinang_use)
    ang2 = np.arccos(cosang_use)

    if ang1 == ang2:
	ang_use = 0.5*(ang1 + ang2)
    else:
        ang1_alt = np.pi - ang1
        if ang1_alt == ang2:
            ang_use = 0.5*(ang1_alt + ang2)
        else:
            ang2_alt = (2*np.pi) - ang2
            if ang1 == ang2_alt:
	        ang_use = 0.5*(ang1 + ang2_alt)
            else:
                ang_use = 0.5*(ang1_alt + ang2_alt)

    l_rad = ang_use + BK
    if l_rad > (2*np.pi):
        l_rad = l_rad - (2*np.pi)

    l_deg = l_rad*180/np.pi
    b_deg = b_rad*180/np.pi

    galcoords = [l_deg, b_deg]
    return galcoords

def get_eq_coords(coordstring):
    stripped = re.sub(r'[^\d.+-]+', '', coordstring)
    #print coordstring
    #print stripped

    rahstring = stripped[0:2]
    ramstring = stripped[2:4]
    rasstring = stripped[4:9]
    decdstring = stripped[9:12]
    decmstring = stripped[12:14]
    decsstring = stripped[14:16]
    pos = False
    if decdstring[0] == '+':
        pos = True

    ra = 15.*(float(rahstring) + (float(ramstring)/60) + (float(rasstring)/3600))
    if pos:
        dec = float(decdstring) + (float(decmstring)/60) + (float(decsstring)/3600)
    else:
        dec = float(decdstring) - (float(decmstring)/60) - (float(decsstring)/3600)

    final_ra = float("{0:.2f}".format(ra))
    final_dec = float("{0:.2f}".format(dec))

#    final_ra = 

    eq_coords = [final_ra, final_dec]
    return eq_coords

def get_totals(subject_vec):
    goodlist = []
    multilist = []
    otherlist = []

    tester_limit = 0

    for entry in subject_vec:
        multifrac = float(entry.classification_dict[u'multi']) / float(entry.num_classifiers)
        goodfrac = float(entry.classification_dict[u'good']) / float(entry.num_classifiers)
        if goodfrac > 0.5:
            entry.majgood = True
	    goodlist.append(entry)
        elif multifrac > 0.5:
            entry.majmult = True
            multilist.append(entry)
        else:
            entry.other = True
            otherlist.append(entry)
        #print entry.zooniverse_id, entry.num_classifiers, entry.classification_dict, multifrac, goodfrac, entry.majmult, entry.majgood, entry.other
        #tester_limit += 1
        #if tester_limit > 5000:
        #    break

    return goodlist, multilist, otherlist


def coord_number_histogram(goodlist, multilist, otherlist, directory):
    #build lists
    goodlat = []
    goodlong = []
    for item in goodlist:
        goodlong.append(item.gal_coords[0])
        goodlat.append(item.gal_coords[1])

    multilat = []
    multilong = []
    for item in multilist:
        multilong.append(item.gal_coords[0])
        multilat.append(item.gal_coords[1])
    
    otherlat = []
    otherlong = []
    for item in otherlist:
        otherlong.append(item.gal_coords[0])
        otherlat.append(item.gal_coords[1])
    
    binvec_lat = np.arange(-90., 91., 5.)

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    plt.hist([goodlat, multilat, otherlat], bins=binvec_lat, stacked=True, label=('Good', 'Multi', 'Other'))
    counts_good,edges_good = np.histogram(goodlat, bins=binvec_lat)
    counts_multi,edges_multi = np.histogram(multilat, bins=binvec_lat)
    counts_other,edges_other = np.histogram(otherlat, bins=binvec_lat)

    f1 = open(directory+'/raw_numbers_binned.dat','w')
    for i in range(0, len(binvec_lat)-1, 1):
        f1.write(str(counts_good[i])+','+str(counts_multi[i])+','+str(counts_other[i])+'\n')
    f1.close()

    totals = np.zeros(len(binvec_lat) - 1)
    rates = np.zeros(len(binvec_lat) - 1)
    rate_errs = np.zeros(len(binvec_lat) - 1)
    for i in range(0,len(binvec_lat)-1,1):
        multi_errs = np.sqrt(float(counts_multi[i]))
        totals[i] += float(counts_good[i] + counts_multi[i] + counts_other[i])
        total_errs = np.sqrt(totals[i])

        if totals[i] > 0.:
            rates[i] += ((float(counts_multi[i]))/float(totals[i]))
            rate_errs[i] = rates[i]*np.sqrt((((multi_errs)**2)/((float(counts_multi[i]))**2)) + ((total_errs**2)/((totals[i])**2)))

    plt.xlabel("Galactic latitude (degrees)")
    plt.xlim([-90., 90.])
    plt.xticks([-90., -60., -30., 0., 30., 60., 90.])
    plt.ylabel("Number of subjects")
    plt.legend(loc = "upper left")

    plt.savefig(directory+"/classification_data_raw_numbers.png")

    print "Number Histogram Plotted"
    return rates, rate_errs

def coord_density_histogram(ratevec, rate_errvec, directory):
    binvec_lat = np.arange(-90., 91., 5.)

    hist_x = []
    hist_y = []
    error_x = []
    error_y = []
 
    for i in range(0,(len(binvec_lat.tolist())-1),1):
        hist_x.append(binvec_lat[i])
        hist_x.append(binvec_lat[i+1])
        error_x.append(0.5*(binvec_lat[i] + binvec_lat[i+1]))

    for val in ratevec:
        hist_y.append(val)
        hist_y.append(val)
        error_y.append(val)

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    plt.figure()
    plt.plot(hist_x, hist_y)
    plt.errorbar(error_x, error_y, yerr=rate_errvec, fmt='b.')
    plt.xlabel("Galactic latitude (degrees)")
    plt.xlim([-90., 90.])
    plt.xticks([-90., -60., -30., 0., 30., 60., 90.])
    plt.ylabel("Contamination rate")

    plt.savefig(directory+"/classification_data_rates.png")
    #print "Histogram saved"
    return "Histogram saved"

def print_to_file(vec, directory):
    f1 = open(directory+'/combined_classification_data_information.csv', 'w')

    f1.write('Zooniverse_ID,WISE_ID,classifiers,good,multi,oval,empty,extended,shift,state\n')

    for entry in vec:
        f1.write(str(entry) + '\n')

    f1.close()

    print "Wrote to file" 

def lat_long_heatmap(goodlist, multilist, otherlist, directory):
    good_lat = []
    good_long = []

    multi_lat = []
    multi_long = []

    other_lat = []
    other_long = []

    total_lat = []
    total_long = []

    for item in goodlist:
        total_lat.append(item.gal_coords[1])
        total_long.append(item.gal_coords[0])

        good_lat.append(item.gal_coords[1])
        good_long.append(item.gal_coords[0])

    for item in multilist:
        total_lat.append(item.gal_coords[1])
        total_long.append(item.gal_coords[0])

        multi_lat.append(item.gal_coords[1])
        multi_long.append(item.gal_coords[0])

    for item in otherlist:
        total_lat.append(item.gal_coords[1])
        total_long.append(item.gal_coords[0])

        other_lat.append(item.gal_coords[1])
        other_long.append(item.gal_coords[0])

    binvec_long = np.arange(0., 361., 5.)
    binvec_lat = np.arange(-90., 91., 5.)

#    long_num = binvec_long.size
#    lat_num = binvec_lat.size

#    long_coords = np.zeros(long_num-1)
#    lat_coords = np.zeros(lat_num-1)

#    for i in range(long_coords.size):
#        long_coords[i] = 0.5*(binvec_long[i] + binvec_long[i+1])

#    for i in range(lat_coords.size):
#        lat_coords[i] = 0.5*(binvec_lat[i] + binvec_lat[i+1])

    plt.figure()
    total_counts, total_xedges, total_yedges, image_total = plt.hist2d(total_long, total_lat, bins=[binvec_long,binvec_lat],norm=mpl.colors.LogNorm())
    #total_array = image_total.get_array()
    #total_verts = image_total.get_offsets()
    plt.xlim([0., 360.])
    plt.ylim([-90., 90.])
    plt.xticks([0., 30., 60., 90., 120., 150., 180., 210., 240., 270., 300., 330., 360.])
    plt.yticks([-90., -60., -30., 0., 30., 60., 90.])
    plt.colorbar()
    plt.savefig(directory+'/total_colormap.pdf')

    plt.figure()
    multi_counts, multi_xedges, multi_yedges, image_multi = plt.hist2d(multi_long, multi_lat, bins=[binvec_long,binvec_lat],norm=mpl.colors.LogNorm())
    #multi_array = image_multi.get_array()
    plt.xlim([0., 360.])
    plt.ylim([-90., 90.])
    plt.xticks([0., 30., 60., 90., 120., 150., 180., 210., 240., 270., 300., 330., 360.])
    plt.yticks([-90., -60., -30., 0., 30., 60., 90.])
    plt.colorbar()
    plt.savefig(directory+'/multi_colormap.pdf')

    print total_counts.shape
    print multi_counts.shape

    ratio_array = np.zeros(multi_counts.shape)

    for i in range(72):
        for j in range(36):
            ratio_array[i][j] = multi_counts[i][j] / total_counts[i][j]

    x2,y2 = np.meshgrid(binvec_long, binvec_lat)

    print (x2.shape, y2.shape)

    print ratio_array.shape
    print np.swapaxes(ratio_array,0,1).shape

    minplane = [-5., -5.]
    maxplane = [5., 5.]

    plt.figure()
    plt.pcolormesh(x2, y2, np.swapaxes(ratio_array,0,1), vmin = 0., vmax = 1.)
    plt.xlim([0., 360.])
    plt.ylim([-90., 90.])
    plt.xticks([0., 30., 60., 90., 120., 150., 180., 210., 240., 270., 300., 330., 360.])
    plt.yticks([-90., -60., -30., 0., 30., 60., 90.])
    plt.colorbar()
    plt.xlabel('Galactic longitude', fontsize=16)
    plt.ylabel('Galactic latitude', fontsize=16)

    plt.plot([0., 360.], minplane, 'w--', linewidth=2)
    plt.plot([0., 360.], maxplane, 'w--', linewidth=2)

    plt.plot([280.4652, 302.8084], [-32.8884, -44.3277], 'wx')

    plt.savefig(directory+'/density_colormap.pdf')

    doneness1 = "Done"

    return doneness1

    #plt.figure()
    #plt.pcolor(
