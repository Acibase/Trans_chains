#####
#	implements the E-model
#
# 	last updated: 10th Oct, 2013

#	NOTE:
#		Works for a single N, at the same location in each sample (representing the plasmid negavive samples)
#		If this is not the case for your data, adapt accordingly
#


import sys
import os
from collections import defaultdict
import operator



psg_dict = defaultdict(list)		# key: patient_id, value: a list of genotypes
clean_psg_dict = defaultdict(list)	# key: patient_id, value: a list of genotypes, see filter_psg_data()
patient_initial = {}			# key: patient_id, value: the ancestral genotype for this patient
patient_stay = defaultdict(list)	# key: patient_id, value list of tuples (ward,att_date,det_date)

reci_donor_dict = defaultdict(list)		# key: recipient_id, value list of tuples (donor_id,dist,day_first_contact)
# sorted on dist, which is the #SNPs between recipient and donor

reci_donor_contact = defaultdict(dict)		# key: recipient_id, 2nd level key: donor_id,
# 2nd level value: day of first contact for this pair

reci_donor_genotype = defaultdict(list)  	# key: recipient_id, value: a list of tuples (donor_id,dist)
# sorted on dist, which is the #SNPs between recipient and donor

reci_donor_final = defaultdict(list)		# key: recipient_id, value list of tuples (donor_id,dist,day_first_contact)
# sorted on dist, and cleaned from donors with first contact day
# later than the first positive day for the recipient


patient_positive_day = {}		# key: patient_id, value: the first positive sample


INDEX_PATIENT = '1'
INDEX_GENOTYPE = 'CTTGGACAAGTGGCTTCGGGTGTGGTACGCT'






def go(psg,stay,positive,days_offset):
    # read the patient_genotype file
    read_psg_data(psg)
    
    # filter it to remove non-informative samples,
    # i.e., if a patient has more than one sample with the same genotype, keep only one
    filter_psg_data()
    
    # pick the initial genotype for each patient:
    #     a) smallest number of SNPs to *index*
    #     b) all other genotypes from this patient can be derived by SNP accumulation
    #     c) cannot be plasmid negative if at least one other sample from this patient is plasmid positive
    #         assumptions: no plasmid gain, progression within host
    pick_initial()
    
    
    # read the patient_stay data
    read_patient_stay(stay)
    
    # based on patient stay ONLY!!!
    # choose potential donors for each patient based on ward overlap
    tot_pot_donors = 0
    patient_keys = patient_initial.keys()
    patient_keys.sort()
    for pk in patient_keys:
        if pk == INDEX_PATIENT:
            continue
        pot_contact_donors = get_contacts(pk,days_offset)
        reci_donor_contact[pk] = pot_contact_donors
        tot_pot_donors = tot_pot_donors + len(pot_contact_donors)
    print "%s patients (excluding index): ON WARD STAY DATA ONLY: Total number of potential donors = %s; %s on average per patient" % (len(psg_dict)-1,tot_pot_donors, tot_pot_donors/(len(patient_initial)-1))
    
    for pk in patient_keys:
        print "%s: %s" % (pk,reci_donor_contact[pk])
    print "=========================================================================="
    
    # based on genomic data ONLY!!!
    # choose potential donors for each patient, for patients with multi samples, use the initial genotype
    tot_pot_donors = 0
    patient_keys = patient_initial.keys()
    patient_keys.sort()
    for pk in patient_keys:
        if pk == INDEX_PATIENT:
            continue
        r_id = pk
        r_geno = patient_initial[r_id]
        pot_geno_donors = find_potential_donors(r_id,r_geno)
        pot_geno_donors_sorted = sorted(pot_geno_donors.iteritems(), key=operator.itemgetter(1))
        reci_donor_genotype[pk] = pot_geno_donors_sorted
        tot_pot_donors = tot_pot_donors + len(pot_geno_donors_sorted)
    print "%s patients (excluding index): ON GENOTYPE DATA ONLY: Total number of potential donors = %s; %s on average per patient" % (len(psg_dict)-1,tot_pot_donors, tot_pot_donors/(len(patient_initial)-1))
    
    for pk in patient_keys:
        print "%s: %s" % (pk,reci_donor_genotype[pk])
    print "=========================================================================="
    
    
    # read patient_positive - the date of the first confirmed outbreak positive sample
    read_positive(positive)
    
    # now, intersect the patient_stay and patient_genotype data!!!
    intersect_stay_genotype()
    
    # now go over the result of intersect_stay_genotype - in reci_donor_dict
    # and for each recipent delete all donors with earliest transmission day > first positive sample
    for reci_id,v in reci_donor_dict.iteritems():
        for d in v:
            donor_id = d[0]
            dist = d[1]
            first_cont = int(d[2])
            if first_cont <= patient_positive_day[reci_id]:
                reci_donor_final[reci_id].append(d)
    
    
    # print out the results
    num_trans = 0
    for k,v in reci_donor_final.iteritems():
        num_trans = num_trans + len(v)
    print "Combining contact and genotype data, checking for first positive day"
    print "Found potential donors for %s patients (out of %s, excl. index)" % (len(reci_donor_final),len(psg_dict)-1)
    print "Total number of donors = %s, %s on average per patient with found donors" % (num_trans,num_trans/len(reci_donor_final))
    for k,v in reci_donor_final.iteritems():
        print "recipient_id = %s, donor data = %s" % (k,v)
    print "=========================================================================="
    
    
    # print out patient_ids for which no donor is found
    for id in clean_psg_dict.keys():
        if id not in reci_donor_final:
            print "No potential donor found for patient_id = %s" % (id)










def read_positive(positive):
    in_han = open(positive,'r')
    for line in in_han:
        pat_id,pos_day = line.strip().split('\t')
        if pat_id not in patient_positive_day:
            patient_positive_day[pat_id] = int(pos_day)
        else:
            print "Duplicate patient_id in %s = %s" % (positive,pat_id)
    in_han.close()
    print "Read the first positive sample data for %s patients" % (len(patient_positive_day))
    print "=========================================================================="




def intersect_stay_genotype():
    for x in reci_donor_genotype.keys():
        reci_id = x
        geno_donor_tuples = reci_donor_genotype[reci_id]
        for gd in geno_donor_tuples:
            dono_id = gd[0]
            dist = gd[1]
            # check if this pair has been in contact
            if dono_id in reci_donor_contact[reci_id]:
                first_contact_day = reci_donor_contact[reci_id][dono_id]
                tuple = (dono_id,dist,first_contact_day)
                reci_donor_dict[reci_id].append(tuple)




def get_contacts(reci_id,days_offset):
    donors = {}         # key: donor's patient_id, value: the earliest day these pair met
    reci_stays = patient_stay[reci_id]  # a list of tuples (ward,att_date,det_date) - all ward stays for this recipient
    for r in reci_stays:
        pk = patient_stay.keys()
        pk.sort()
        for dono_id in pk:
            if dono_id == reci_id:
                continue
            else:
                dono_stays = patient_stay[dono_id]	# a list of tuples (ward,att_date,det_date) - all ward stays for this donor
                for d in dono_stays:
                    if r[0] != d[0]:    # these two stays for these two patients were on diff wards
                        continue
                    if r[1] > d[2]+days_offset:     # recipient attached after donor detached (+OFFSET!!!)
                        continue
                    if r[2] < d[1]:	# recipient detached before donor attached
                        continue
                    # the recipent and donor overlaped at this ward in this time inetrval
                    # establish day of first contact for this pair
                    if r[1] > d[1]:
                        first_contact_day = r[1]
                    else:
                        first_contact_day = d[1]
                    if dono_id not in donors:
                        donors[dono_id] = first_contact_day
                    else:
                        day_recorded = donors[dono_id]
                        if first_contact_day < day_recorded:
                            donors[dono_id] = first_contact_day
    return donors



def read_patient_stay(stay):
    in_han = open(stay,'r')
    for line in in_han:
        data = line.strip().split('\t')
        pat_id = data[0].strip()
        for x in xrange(1,len(data),3):
            ward = data[x].strip()
            att_date = int(data[x+1])
            det_date = int(data[x+2])
            if att_date > 1000:
                print "Out of scope att_date: pat_id = %s, ward = %s, att_date = %s, det_date = %s" % (pat_id,ward,att_date,det_date)
                raise SystemExit
            if det_date > 1000:
                print "Out of scope det_date: pat_id = %s, ward = %s, att_date = %s, det_date = %s" % (pat_id,ward,att_date,det_date)
                raise SystemExit
            tuple = (ward,att_date,det_date)
            patient_stay[pat_id].append(tuple)
    
    tot_ward_stays = 0
    for k,v in patient_stay.iteritems():
        tot_ward_stays = tot_ward_stays + len(v)
    
    print "=========================================================================="
    print "Patient stay data read: %s patients, %s ward stays" % (len(patient_stay),tot_ward_stays)



def find_potential_donors(r_id, r_genotype):
    donors = {}				# key: donor_id, value: the smallest distance
    for d_id,list_d_genotypes in clean_psg_dict.iteritems():
        if d_id == r_id:		# that's the same patient
            continue
        for d_genotype in list_d_genotypes:
            dist = compute_distance(r_genotype,d_genotype)
            if dist < 0:
                continue
            else:
                if d_id not in donors:
                    donors[d_id] = dist
                else:				# make sure that we get the closest genotype from each potential donor
                    dist_there = donors[d_id]
                    if dist < dist_there:
                        donors[d_id] = dist
    return donors


def compute_distance(r,d):
    # r - recepient
    # d - donor
    # use the index case genotype
    ref = INDEX_GENOTYPE
    if len(r) != len(d):
        print "Length mismatch in compute_distance, r = %s, d = %s" % (r,d)
        raise SystemExit
    dist = 0
    for x in xrange(0,len(ref)):
        if r[x] != d[x]:    	        # a SNP between the two
            if d[x] == 'N':          	# this sample lost the plasmid, assumption no plasmid gain - cannot be donor
                dist = -1
                break
            elif r[x] == ref[x]:	# recepient same as index, donor acquired a SNP - a flip back, not allowed
                dist = -1
                break
            else:			# a SNP in recepient compared to the donor
                dist += 1
    return int(dist)


def pick_initial():
    num_single = 0
    num_multi = 0
    for k,v in clean_psg_dict.iteritems():
        if len(v) == 1:				# single genotype observed in this patient, it is the initial one
            patient_initial[k] = v[0]
        else:
            chosen_genotype = pick_one_genotype(k,v)
            patient_initial[k] = chosen_genotype
    print "Choosing the initial genotype for each patient done"
#    for k,v in patient_initial.iteritems():
#        print "Patient %s: %s" % (k,v)


def pick_one_genotype(k,genotypes):
    dist_to_ref = {}	# key: genotype, value: distance (num of SNPs)
    for g in genotypes:
        dist = int(get_distance(g))
        dist_to_ref[g] = dist
    sorted_dist = sorted(dist_to_ref.iteritems(), key=operator.itemgetter(1))
    
    # pick the one closest to the index genotype but not plasmid negative, unless all are plasmid negative
    chosen_initial = pick_closest_noN(sorted_dist)
    if chosen_initial == 'None':
        print "It was not possible to choose the initial genotypes for Patient %s" % (k)
        raise SystemExit
    
    #    print "Patient %s:" % (k)
    #    print sorted_dist
    #    print "chosen initial genotype = %s" % (chosen_initial)
    
    # now check if the remaining genotypes can be derived from the initial, without any flip-backs
    all_fine = True
    for x in xrange(1,len(sorted_dist)):	# skipping the first one, that's the pot_initial
        ok = can_derive(chosen_initial,sorted_dist[x][0])
        if not ok:
            all_fine = False
            print "It was not possible to derive all other genotypes from the pot_initial genotype for Patient %s" % (k)
            raise SystemExit
    return chosen_initial


def pick_closest_noN(sorted_dist):
    chosen_genotype = 'None'
    for x in xrange(0,len(sorted_dist)):
        pot_initial = sorted_dist[x][0] 	# i.e., start with the one closest to the index genotype
        if pot_initial.find('N') == -1:		# and is plasmid positive
            chosen_genotype = pot_initial
            break
        else:					# is plasmid negative
            # check if all genotypes from this patient are also plasmid negative
            all_N = True
            for y in xrange(x+1,len(sorted_dist)):
                if sorted_dist[y][0].find('N') == -1:
                    all_N = False
                    break
            if all_N:				# the pot_initial is plasmid negative, but so are all the others
                chosen_genotype = pot_initial
                break
    # else, pick the next one (in terms of distance to index genotype)
    return chosen_genotype


def can_derive(sour, dest):
    can = True
    # use the index case genotype
    ref = INDEX_GENOTYPE
    for x in xrange(0,len(ref)):
        if sour[x] != dest[x]:		# a SNP between the two, source must be same as ref at this location
            if dest[x] == 'N':		# this sample lost the plasmid, ok
                continue
            if sour[x] != ref[x]:
                can = False
                print "x = %s, sour[x] = %s, dest[x] = %s, ref[x] = %s" % (x,sour[x],dest[x],ref[x])
                break
    return can


def get_distance(g):
    dist = 0
    # use the index case genotype
    ref = INDEX_GENOTYPE
    if len(ref) != len(g):
        print "Length mismatch between index genotype = %s and g = %s" % (ref,g)
        raise SystemExit
    for x in xrange(0,len(ref)):
        if ref[x] != g[x]:
            dist += 1
    return dist


def filter_psg_data():
    num_dupl_genotypes = 0
    for k,v in psg_dict.iteritems():
        if len(v) == 1:
            clean_psg_dict[k] = v
        else:
            all_genotypes = v
            for g in all_genotypes:
                if k not in clean_psg_dict:		# seeing this patient for the first time
                    clean_psg_dict[k].append(g)
                else:
                    list_clean_genotypes = clean_psg_dict[k]
                    if g in list_clean_genotypes:
                        #                        print "Patient %s, Duplicate genotype = %s" % (k,g)
                        num_dupl_genotypes += 1
                    else:
                        clean_psg_dict[k].append(g)
    print "Patient genotype filtering done, found a total of %s duplicate genotypes" % (num_dupl_genotypes)
    clean_samples = 0
    for k,v in clean_psg_dict.iteritems():
        clean_samples = clean_samples + len(v)
    print "Found %s patients, total of %s patient-unique samples (including index)" % (len(psg_dict),clean_samples)
    sys.stdout.flush()


def read_psg_data(psg):
    in_han = open(psg,'r')
    for line in in_han:
        data = line.strip().split('\t')
        pat_id = data[0].strip()
        geno_end = len(data)
        genotype = ''
        for x in xrange(1,geno_end):
            genotype = genotype + "%s" % (data[x])
        psg_dict[pat_id].append(genotype)
    
    tot_samples = 0
    for k,v in psg_dict.iteritems():
        tot_samples = tot_samples + len(v)
    print "Patient Genotype data read"
    print "Found %s patients, total of %s samples (including index)" % (len(psg_dict),tot_samples)
    sys.stdout.flush()










if __name__ == '__main__':
    if (len(sys.argv) != 5) :
        print '''********************\n
            Usage suggestion: time python trans_chains.py patient_genotype.txt patient_stay.txt patient_positive.txt days_offset-as an integer\n
            python \n
            *********************'''
    else:
        go(str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[3]),int(sys.argv[4]))
