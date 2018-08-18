__author__ = 'Michael'

import random
import os, sys
import dendropy
import numpy as np
import platform
import datetime as dt
import multiprocessing

if platform.system()=='Darwin':
    db='/Users/michaelnute/Dropbox/KRA Primate Project/16S'
    work='/Users/michaelnute/Dropbox/Grad School/Phylogenetics/work/kra-primate-project/kra-primate'
elif platform.system()=='Windows':
    db = 'C:/Users/Michael/Dropbox/KRA Primate Project/16S'
    work='C:/Users/Michael/Dropbox/Grad School/Phylogenetics/work/kra-primate-project/kra-primate'
elif platform.system()=='Linux':
    work='/projects/tallis/nute/work/metagenomics/kra-primate'
    db=work

npdata=os.path.join(work,'npdata')
# tnew=dendropy.Tree.get(file=open(os.path.join(npdata,'observation_id_subtree.tre'),'r'),schema="newick",preserve_underscores=True)
tnew=dendropy.Tree.get(file=open(os.path.join(npdata,'rep_set_pasta.tre'),'r'),schema="newick",preserve_underscores=True)

def make_D_matrix_parallel():
    obs_ids=open(os.path.join(npdata,'obs-ids.txt'),'r')
    taxa={}
    ct=0
    for i in obs_ids:
        a=i.strip()
        taxa[a]=ct
        ct+=1
    lt=len(taxa)
    obs_ids.close()
    print (lt)
    # dmat=np.zeros((lt,lt),dtype=np.float64)

    qs = []
    for i in range(16):
        qs.append(multiprocessing.Queue())

    args=[]
    ct = 0
    for i in tnew.preorder_node_iter():
        if i.edge_length != None:
            # tic=dt.datetime.now()
            l=i.edge_length
            mytaxa=[]
            for a in i.leaf_iter():
                mytaxa.append(a.taxon.label)
            args.append((l,mytaxa))
            ind = ct % 16
            qs[ind].put((   l,mytaxa))

            # toc=dt.datetime.now()-tic
            # print "\t%s" % str(toc)
            print (ct)
            ct +=1
    print ("done making args list")

    p00=multiprocessing.Process(target=worker,args=(qs[0],'q0',taxa.copy()))
    p00.start()
    p01=multiprocessing.Process(target=worker,args=(qs[1],'q1',taxa.copy()))
    p01.start()
    p02=multiprocessing.Process(target=worker,args=(qs[2],'q2',taxa.copy()))
    p02.start()
    p03=multiprocessing.Process(target=worker,args=(qs[3],'q3',taxa.copy()))
    p03.start()
    p04=multiprocessing.Process(target=worker,args=(qs[4],'q4',taxa.copy()))
    p04.start()
    p05=multiprocessing.Process(target=worker,args=(qs[5],'q5',taxa.copy()))
    p05.start()
    p06=multiprocessing.Process(target=worker,args=(qs[6],'q6',taxa.copy()))
    p06.start()
    p07=multiprocessing.Process(target=worker,args=(qs[7],'q7',taxa.copy()))
    p07.start()
    p08=multiprocessing.Process(target=worker,args=(qs[8],'q8',taxa.copy()))
    p08.start()
    p09=multiprocessing.Process(target=worker,args=(qs[9],'q9',taxa.copy()))
    p09.start()
    p10=multiprocessing.Process(target=worker,args=(qs[10],'q10',taxa.copy()))
    p10.start()
    p11=multiprocessing.Process(target=worker,args=(qs[11],'q11',taxa.copy()))
    p11.start()
    p12=multiprocessing.Process(target=worker,args=(qs[12],'q12',taxa.copy()))
    p12.start()
    p13=multiprocessing.Process(target=worker,args=(qs[13],'q13',taxa.copy()))
    p13.start()
    p14=multiprocessing.Process(target=worker,args=(qs[14],'q14',taxa.copy()))
    p14.start()
    p15=multiprocessing.Process(target=worker,args=(qs[15],'q15',taxa.copy()))
    p15.start()


def worker(q,qid,taxadict):
    lt=len(taxadict)
    dmat=np.zeros((lt,lt),dtype=np.float64)
    while q.empty()==False:
        # try:
        item=q.get()
        l=item[0]
        tinds=[taxadict[i] for i in item[1]]
        for i in range(len(tinds)):
            ind1=tinds[i]
            dmat[ind1,ind1]+=l
            for j in range(i+1,len(tinds)):
                ind2=tinds[j]
                dmat[ind1,ind2]+=l
                dmat[ind2,ind1]+=l
    np.savetxt(os.path.join(work,'dmatrices','dMatrix-' + qid + '.csv'),dmat,delimiter=',')
    print ("done with %s" % qid)

def assemble_d_matrices():
    nt=6862
    dmat=np.zeros((nt,nt),dtype=np.float64)
    for i in range(16):
        name=os.path.join(work,'dmatrices','dMatrix-')
        name = name + 'q%s.csv' % i
        dmatnew=np.loadtxt(name,dtype=np.float64,delimiter=',')
        dmat+=dmatnew
        del dmatnew
        print ("done %s" % i)
    np.savetxt(os.path.join(work,'dmatrices','dMatrix-full.csv'),dmat,delimiter=',')
    outf=open(os.path.join(npdata,'dmatrixfull.bin'),'wb')
    dmat.tofile(outf)
    outf.close()
    print ("Wrote full dmatrix in csv form to:\n\t %s" % os.path.join(work,'dmatrices','dMatrix-full.csv'))
    print ("and in binary form to:\n\t%s" % os.path.join(npdata,'dmatrixfull.bin'))
    return dmat

def get_taxa_indices():
    obs_ids=open(os.path.join(npdata,'obs-ids.txt'),'r')
    taxa={}
    ct=0
    for i in obs_ids:
        a=i.strip()
        taxa[a]=ct
        ct+=1
    lt=len(taxa)
    obs_ids.close()
    print (lt)
    return taxa

def make_D_matrix():

    # return tnew
    import datetime
    obs_ids=open(os.path.join(npdata,'obs-ids.txt'),'r')
    taxa=[]
    for i in obs_ids:
        a=i.strip()
        taxa.append(a)
    lt=len(taxa)
    print (lt)
    dmat=np.zeros((lt,lt),dtype=np.float64)

    tic=datetime.datetime.now()
    for i in range(lt):
        # if i % 100==0:
        # print "i = %s" % i
        for j in range(i,lt):
            if j % 1000 ==0:
                print ("\tj = %s" % j)
            l=tnew.mrca(taxon_labels=set([taxa[i],taxa[j]])).distance_from_root()
            dmat[i,j]=l
            dmat[j,i]=l
        toc=datetime.datetime.now()

        print ("i = %s: %s" % (i,str(toc-tic)))

    np.savetxt(os.path.join(npdata,'d-matrix.csv'),dmat,delimiter=',')
    outf = open(os.path.join(npdata,"dmatrixfull.bin"),'wb')
    dmat.tofile(outf)
    outf.close()
    print ("done making d-matrix")

    dL, dQ = np.linalg.eig(dmat)
    f=open(os.path.join(npdata,'dmatrixEigenQ.bin'),'wb')
    dQ.tofile(f)
    f.close()
    f1=open(os.path.join(npdata,'dmatrixEigenVals.bin'),'wb')
    dL.tofile(f1)
    f1.close()

def make_DL_DQ(dmat=None):
    if dmat==None:
        df = open(os.path.join(npdata,"dmatrixfull.bin"),'rb')
        dmat = np.fromfile(df,dtype=np.float64)
        df.close()

    dL, dQ = np.linalg.eig(dmat)
    f=open(os.path.join(npdata,'dmatrixEigenQ.bin'),'wb')
    dQ.tofile(f)
    f.close()
    f1=open(os.path.join(npdata,'dmatrixEigenVals.bin'),'wb')
    dL.tofile(f1)
    f1.close()
    print ("write dL (eigen vals) and dQ (eigen vecs) to\n\t%s\nand\n\t%s\nrespectively" % (os.path.join(npdata,'dmatrixEigenVals.bin'),os.path.join(npdata,'dmatrixEigenQ.bin')))


def test_pair(a=None):
    if a==None:
        a=random.sample(taxa.keys(),2)
    print (a)
    print (tnew.mrca(taxon_labels=a).distance_from_root())
    # print (dmat[taxa[a[0]],taxa[a[1]]])

taxa=get_taxa_indices()

if __name__=='__main__':
    if '--makedmat' in sys.argv:
        # make_D_matrix()
        # make_D_matrix_parallel()
        dm=assemble_d_matrices()
        make_DL_DQ(dm)
    else:
        print ("no argument specified") 
    # test_leaf_nodes()
    # test_pair()