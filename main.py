__author__ = 'Michael'
import os
import numpy as np
# import scipy as sp
# import biom
import dendropy
import platform, multiprocessing, datetime
# from guppy import hpy

if platform.system()=='Darwin':
    db='/Users/michaelnute/Dropbox/KRA Primate Project/16S'
    work='/Users/michaelnute/Dropbox/Grad School/Phylogenetics/work/kra-primate-project/kra-primate'
elif platform.system()=='Windows':
    db = 'C:/Users/Michael/Dropbox/KRA Primate Project/16S'
    work='C:/Users/Michael/Dropbox/Grad School/Phylogenetics/work/kra-primate-project/kra-primate'
elif platform.system()=='Linux':
    work='/projects/tallis/nute/work/metagenomics/kra-primate'
    db=work
# tree_loc = os.path.join(db,'rep_set.tre')
tree_loc = os.path.join(work,'npdata','rep_set_pasta.tre')
biom_loc = os.path.join(work,'otu_table_w_tax_open_ref_r1r2r3_nlf_primates_nomac_noav_noprop_novar.biom')

npdata=os.path.join(work,'npdata')
# xtable_loc = os.path.join(work,'primates-raw-xtable.csv')
# xtable_loc = os.path.join(work,'primates-raw-xtable-5vars.csv')
# xtable_loc = os.path.join(work,'primates-raw-xtable-8vars.csv')
xtable_loc = os.path.join(npdata,'x-data.csv')
ytable_loc = os.path.join(npdata,'y-data.csv')
dQ_loc = os.path.join(npdata,'dmatrixEigenQ.bin')
dL_loc = os.path.join(npdata,'dmatrixEigenVals.bin')
dFull_loc = os.path.join(npdata,'dmatrixfull.bin')

def initurs():
    urs=UnifracRegressionSolver(ytable_loc,xtable_loc,tree_loc,None,dQ_loc,dL_loc,dFull_loc)
    urs.initialize_all()
    return urs



class UnifracRegressionSolver():
    data=None
    y=None
    yhat=None
    ehat=None
    x=None
    Beta=None
    tree=None
    dMatrix=None
    lda=1
    stat_T=0
    stat_P=0
    stat_N=0
    jackknifes={}
    numprocs=12

    def __init__(self,yfile=None, xfile=None, treefile=None, outputfile=None, dEigenQ=None, dEigenL=None, dfile=None, obs_ids_file=None,sample_ids_file=None):

        self.yfile=yfile
        self.xfile=xfile
        self.treefile=treefile
        self.dfile=dfile
        self.outfile=outputfile
        self.obs_ids_f=obs_ids_file
        self.samp_ids_f=sample_ids_file
        self.dq_loc = dEigenQ
        self.dl_loc = dEigenL

    def initialize_all(self):
        self.init_make_x()
        self.init_make_y()
        self.init_make_beta_matrix()
        self.init_make_d_matrix(self.dq_loc,self.dl_loc)

    def init_make_x(self, xfile=None):
        '''
        Reads a CSV file with the x-matrix and converts it to a numpy array, then re-reads it to read the labels along
        the left-hand column into the variable 'xlabs'
        :param xfile: path to to the x-matrix file
        :return: None
        '''
        if xfile is not None:
            self.xfile=xfile

        self.x = np.loadtxt(self.xfile, np.float64, delimiter=',')
        self.stat_N = self.x.shape[1]
        self.stat_P = self.x.shape[0]
        print 'X matrix is (%s x %s) and has been loaded from:' % (str(self.stat_P),str(self.stat_N))
        print self.xfile + '\n'


    def init_make_y(self,yfile=None):
        if yfile<>None:
            self.yfile=yfile
        self.y=np.loadtxt(self.yfile,dtype=np.float64,delimiter=',')
        self.stat_T = self.y.shape[0]
        print 'Y matrix is (%s x %s) and has been loaded from:' % (str(self.y.shape[0]),str(self.y.shape[1]))
        print self.yfile + '\n'


    def init_make_d_matrix(self,dEigenQ=None,dEigenVals=None):
        '''
        imports the D matrix if the Q and V have not already been calculated.
        :return:
        '''
        self.dMatrix=np.fromfile(open(self.dfile,'rb'),dtype=np.float64).reshape(self.stat_T,self.stat_T)
        if dEigenQ==None or dEigenVals==None:
            print '''
            The Eigendecomposition for the D-matrix have not been provided, so you will need
                to execute the method 'self.init_run_Dmatrix_eigendecomp()' before proceeding
                      '''
        else:
            self.dQ=np.fromfile(open(dEigenQ,'rb'),dtype=np.float64).reshape(self.stat_T,self.stat_T)
            self.dLamda=np.fromfile(open(dEigenVals,'rb'),dtype=np.float64).reshape(self.stat_T,1)
            print 'D-Matrix Eigendecomposition successfully loaded'

    def init_run_Dmatrix_eigendecomp(self):
        self.dQ,self.dLamda = np.linalg.eig(self.dMatrix)

    def init_make_beta_matrix(self):
        '''
        initializes the matrix of beta values
        :return:
        '''
        self.Beta=np.zeros((self.stat_T,self.stat_P),np.float64)

    def predict(self,ytmp,xtmp,btmp):
        '''
        puts B*x into the value of yhat. clips it at a min of zero (ACTUALLY NO)
        :return:
        '''
        # np.clip(np.dot(self.Beta,self.x),0,np.inf,self.yhat)

        #change this part when we override this method
        yhat = np.dot(btmp,xtmp)
        ehattmp = ytmp-yhat

        #return squared unifracs:
        my_N=ytmp.shape[1]
        EvansMatsen=np.zeros(my_N,dtype=np.float64)

        for i in range(my_N):
            EvansMatsen[i]=np.dot(np.dot(ehattmp[:,i].T,self.dMatrix),ehattmp[:,i].T)
        return ehattmp, EvansMatsen

    def tofd(self,ld=1):
        self.train_on_full_data(ld)

    def train_on_full_data(self,lamda=None):
        if lamda<>None:
            self.lda=lamda
        print '%s: training the model on the full data with lambda = %s.' % (str(datetime.datetime.now()),str(self.lda))
        self.train(self.y,self.x,self.lda,self.dQ,self.dLamda,self.stat_T,self.stat_P)
        print '%s: Done training.' % str(datetime.datetime.now())

    def train(self,*args):
        self.Beta=train_linear(*args)

    def run_cv_for_lamda(self,lamda):
        print '%s: Running leave one out error for lambda = %s' % (str(datetime.datetime.now()),str(lamda))
        residuals=np.ones(self.stat_N,dtype=np.float64)*(-1)
        # refs=(self.stat_N,self.x,self.y,self.dQ,self.dLamda,self.dMatrix,self.stat_T,self.stat_P)


        dc=[]
        ppipes=[]
        procs=[]
        p_results=[]
        for i in range(self.numprocs):
            dc.append(self.get_data_copies())
            (ps,pr)=multiprocessing.Pipe()
            p1=multiprocessing.Process(target=leave_one_out_thread_worker,args=(range(self.stat_N)[i::self.numprocs],lamda,dc[i],ps,i))
            p1.start()
            procs.append(p1)
            ppipes.append((ps,pr))

            # datacopy1=self.get_data_copies()
            # p1_s, p1_r=multiprocessing.Pipe()
            # p1 = multiprocessing.Process(target=leave_one_out_thread_worker,args=(range(self.stat_N)[::4],lamda,datacopy1,p1_s,1))
            # p1.start()

        # datacopy2=self.get_data_copies()
        # p2_s, p2_r=multiprocessing.Pipe()
        # p2 = multiprocessing.Process(target=leave_one_out_thread_worker,args=(range(self.stat_N)[1::4],lamda,datacopy2,p2_s,2))
        # p2.start()
        #
        # datacopy3=self.get_data_copies()
        # p3_s, p3_r=multiprocessing.Pipe()
        # p3 = multiprocessing.Process(target=leave_one_out_thread_worker,args=(range(self.stat_N)[2::4],lamda,datacopy3,p3_s,3))
        # p3.start()
        #
        # datacopy4=self.get_data_copies()
        # p4_s, p4_r=multiprocessing.Pipe()
        # p4 = multiprocessing.Process(target=leave_one_out_thread_worker,args=(range(self.stat_N)[3::4],lamda,datacopy4,p4_s,4))
        # p4.start()

        for i in range(self.numprocs):
            p_res=ppipes[i][1].recv()
            p_results.append(p_res)

        # p1_results = p1_r.recv()
        # p2_results = p2_r.recv()
        # p3_results = p3_r.recv()
        # p4_results = p4_r.recv()

        for pres in p_results:
            for i in pres.keys():
                residuals[i]=pres[i]

        # for i in p1_results.keys():
        #     residuals[i]=p1_results[i]
        # for i in p2_results.keys():
        #     residuals[i]=p2_results[i]
        # for i in p3_results.keys():
        #     residuals[i]=p3_results[i]
        # for i in p4_results.keys():
        #     residuals[i]=p4_results[i]

        print sum(residuals>=0)
        self.jackknifes[lamda]=residuals
        print "%s: done running the leave one out error for lambda, residuals: %s" % (str(datetime.datetime.now()),str(sum(residuals)))

    def run_cv_for_lambda_grid(self,values=None):
        if values==None:
            values=[0,.0001,.001,.01,.05,.075,.1,.25,.5,.75,1,5,10]
        for i in values:
            self.run_cv_for_lamda(i)
            # self.save_beta_values(os.path.join(npdata,'beta_pasta'),'beta_pasta_lambda-%s.txt' % i)

    def get_data_copies(self):
        refs={
            'x_cp':self.x.copy(),
            'y_cp':self.y.copy(),
            'N':self.stat_N,
            'Tx':self.stat_T,
            'P':self.stat_P,
            'dQ':self.dQ.copy(),
            'dL':self.dLamda.copy(),
            'dMat':self.dMatrix.copy()
        }
        return refs

    def make_output_files(self,output_folder=None,outputfile=None):
        res=None
        for i in self.jackknifes.keys():
            a=np.hstack((np.array([i],dtype=np.float64),self.jackknifes[i]))
            if res==None:
                res=a
            else:
                res=np.vstack((res,a))
        if output_folder==None:
            output_folder='./'
        if outputfile==None:
            output_file='residuals.csv'
        outpath=os.path.join(output_folder,outputfile)
        np.savetxt(outpath,res,delimiter=',')

    def save_beta_values(self,output_folder=None,outputfile=None):
        if output_folder==None:
            output_folder='./'
        if outputfile==None:
            outputfile='beta_outputs.csv'
        outpath=os.path.join(output_folder,outputfile)
        np.savetxt(outpath,self.Beta,delimiter=',')


def leave_one_out_thread_worker(index_todo,lamda,refs,p_sender,thread_id):
    # i=pool_args[0]
    # lamda=pool_args[1]
    # refs=pool_args[2]

    stat_N=refs['N']
    x=refs['x_cp']
    y=refs['y_cp']
    myDQ=refs['dQ']
    myDL=refs['dL']
    dMatrix=refs['dMat']
    Tx=refs['Tx']
    P=refs['P']
    results={}

    for i in index_todo:
        inds=np.arange(stat_N)[np.arange(stat_N)<>i]
        tempx = x[:,inds].copy()
        tempy = y[:,inds].copy()
        # myDfull=self.dMatrix.copy()

        hold_x=x[:,i].copy()
        hold_y=y[:,i].copy()
        myb = train_linear(tempy,tempx,lamda,myDQ,myDL,Tx,P)

        res, emloss = predict_linear(hold_y,hold_x,myb,dMatrix)
        if (i-thread_id) % 100 == 0:
            print 'L-O-O sample %s -- by thread %s -- lambda %s' % (str(i),str(thread_id),str(lamda))
            # q.put((i,emloss))
        results[i]=np.asscalar(emloss)

    p_sender.send(results)


def predict_linear(ytmp,xtmp,btmp,D):
    yhat = np.dot(btmp,xtmp)
    ehattmp = ytmp-yhat

    #return squared unifracs:
    if len(ytmp.shape)>1:
        my_N=ytmp.shape[1]
    else:
        my_N=1
    EvansMatsen=np.zeros(my_N,dtype=np.float64)

    if my_N>1:
        for i in range(my_N):
            EvansMatsen[i]=np.dot(np.dot(ehattmp[:,i].T,D),ehattmp[:,i].reshape((ytmp.shape[0],1)))
    else:
        EvansMatsen[0]=np.dot(np.dot(ehattmp.T,D),ehattmp)
    return ehattmp, EvansMatsen

def numpy_kronecker_mvm(A,B,c,debug=False):
    if debug==True:
        print 'c shape: %s' % str(c.shape)
        print 'A shape: %s' % str(A.shape)
        print 'B shape: %s' % str(B.shape)
    Aht = A.shape[0]
    Awd = A.shape[1]
    Bht = B.shape[0]
    Bwd = B.shape[1]
    bcprime=np.transpose(np.dot(B,np.transpose(np.reshape(c,(Awd,Bwd)))))
    out=np.reshape(np.dot(A,bcprime),(Awd*Bwd,1))
    return out

def train_linear(ytmp,xtmp,lamda,myDQ,myDL,Tx,P,debug=False):
    '''
    Trains the model by fitting the parameters (betas). This can be overridden for a different type of model.
    This function will assume that ytmp, xtmp and btmp are all copies, and it will modify btmp in place.
    :param ytmp:
    :param xtmp:
    :param btmp:
    :param leave_out:
    :return:
    '''
    if debug==True:
        print "Training the linear Unifrac Regression Model"

    u,s,v = np.linalg.svd(xtmp,False)
    my_stat_N = ytmp.shape[1]
    tmp=np.zeros((Tx*P,1),dtype=np.float64)

    #Step 1: Get the last vector y \otimes x
    if debug==True:
        print 'Step 1: Get the last vector y \otimes x'
    btmp=np.zeros((Tx*P,1),dtype=np.float64)
    for i in range(my_stat_N):
        btmp+=np.kron(ytmp[:,i].reshape((Tx,1)),xtmp[:,i].reshape(P,1))
    if debug==True:
        print 'btmp: %s entries nonzero' % str(sum(sum(btmp>0)))

    #Step 2: First KMVM
    if debug==True:
        print 'Step 2: First KMVM'
    tmp=numpy_kronecker_mvm(myDQ.T,u.T,btmp)
    if debug==True:
        print 'tmp shape after step 2: %s' % str(tmp.shape)

    #Step 3: Mutliplying by the diagonal matrix
    if debug==True:
        print 'Step 3: Mutliplying by the diagonal matrix'
    btmp=np.multiply(np.multiply(1/(lamda+np.kron(myDL,np.multiply(s,s))),np.kron(myDL,np.ones(P))).reshape((P*Tx,1)),tmp)

    #Step 4: Second KMVM
    if debug==True:
        print 'Step 4: Second KMVM'
    tmp=numpy_kronecker_mvm(myDQ,u,btmp)

    btmp=tmp.reshape((Tx,P))

    # del myDQ
    # del myDL
    del u, v
    del tmp
    return btmp



# urs.tofd()
print '\n\n'
if __name__=='__main__':
    urs = initurs()
    urs.train_on_full_data(lamda=0.05)
    urs.save_beta_values(npdata,'beta_pastatree_lambda_0.05.txt')
    # urs.run_cv_for_lambda_grid()
    # urs.run_cv_for_lambda_grid([.05,.075])


