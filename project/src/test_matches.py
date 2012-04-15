# To change this template, choose Tools | Templates
# and open the template in the editor.

import numpy as np
import numpy.random as rnd
import numpy.linalg as alg
import geom
import computePose
import solveHomography as sh
import os
import pickle

if __name__ == '__main__':

    set = 'oakland-set'
    query = '2011-10-28_11-51-29_558'
    set2 = 'oakland-unset'
    query2 = query

    ######################

    matches_file = os.path.join('/media/DATAPART2/ah/pose_runs',set,'matches_'+query+'.pkl')
    matches_data = open(matches_file,'rb')
    matches = pickle.load(matches_data)
    matches_data.close()

    matches_file = os.path.join('/media/DATAPART2/ah/pose_runs',set2,'matches_'+query2+'.pkl')
    matches_data = open(matches_file,'rb')
    matches2 = pickle.load(matches_data)
    matches_data.close()


    wRq, wRd, qYaw, nYaw = matches['constants']
    prm1, prm2, imask = matches['iprm'], matches2['iprm'], matches2['imask']

    qray, dray, ddep = matches['qray'][imask,:], matches['dray'][imask,:], matches['ddep'][imask,:]
#    nYaw = prm2[3]
    constants = (wRq,wRd,qYaw,nYaw)
    prm = sh.lsqH_dt(prm1,qray,dray,ddep,constants)
    errs = sh.errH_dt(prm,qray,dray,ddep,constants)
    errs2 = sh.errH_dtn(prm2,qray,dray,ddep,constants)

    print prm[3]*prm[:3]

    print errs
    print np.sum(errs<.01)

    ipool = np.nonzero(imask)[0]
    iq, id, idep = qray, dray, ddep
    print ipool
    
    qray, dray, ddep, qidx, weights = matches['qray'], matches['dray'], matches['ddep'], matches['qidx'], matches['weight']
    nmat, numq = matches['nmat'], matches['numq']
    prm_good = prm
    
    nvalid = 0
    for i in range(1):
        i1 = rnd.randint(0,56)
        i2 = rnd.randint(0,55)
        i2 = i2+1 if i2>=i1 else i2
        i1, i2 = ipool[i1], ipool[i2]
        i1, i2 = 0, 16
        q1,d1,dep1 = qray[i1,:], dray[i1,:], ddep[i1]
        q2,d2,dep2 = qray[i2,:], dray[i2,:], ddep[i2]
        q,d,dep = np.array([q1,q2]), np.array([d1,d2]), np.array([dep1,dep2])
        p1, valid = sh.compH_dt(q,d,dep,constants)
        tmp1 = np.array([0,0,0,1])
        tmp2 = np.array([0,0,0,0,1])
        p1 = sh.lsqH_dt(tmp1,q,d,dep,constants)
        p2, valid = sh.compH_dtn(q,d,dep,constants)
        p2 = sh.lsqH_dtn(p2,q,d,dep,constants)
        

#        print prm_good
#        print prm2
#        print prm1
        errs1 = sh.errH_dt(p1,iq,id,idep,constants)
        errs2 = sh.errH_dtn(p2,iq,id,idep,constants)
        errs_test = sh.errH_dt(p1,q,d,dep,constants)

        print errs1
        print errs2
        print errs_test
        print nYaw
        print p1
        print p2
        print prm_good
        imask, numi, iconf = sh.getInliers(errs,weights,.01,qidx,numq,nmat)
        nvalid += valid



