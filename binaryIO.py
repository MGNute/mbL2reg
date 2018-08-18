__author__ = 'michaelnute'

testfile='testfile.foo'

import os, platform
import numpy as np
import datetime as dt

if platform.system()=='Darwin':
    db='/Users/michaelnute/Dropbox/KRA Primate Project/16S'
    work='/Users/michaelnute/Dropbox/Grad School/Phylogenetics/work/kra-primate-project/kra-primate'
elif platform.system()=='Windows':
    db = 'C:/Users/Michael/Dropbox/KRA Primate Project/16S'
    work='C:/Users/Michael/Dropbox/Grad School/Phylogenetics/work/kra-primate-project/kra-primate'
    csvbackup='C:\\Users\\Michael\\Grad School Stuff\\Research\\Phylogenetics\\Metagenomics\\kra-primate-csv'
elif platform.system()=='Linux':
    work='/projects/tallis/nute/work/metagenomics/kra-primate'
    db=work

# foo=np.array([[2.87,4.23,5235.23,12.4,123.4],[ 1,2,3,4,5]],dtype=np.float64)
#
# #foo=2.87
#
# fout= open(os.path.join(workfolder,'project',testfile),'wb')
# foo.tofile(fout) #stores it in row major order (across, then down)
# fout.close()

npdata=os.path.join(work,'npdata')

def convert_numpy_csv_to_binary(inpath,outpath,rowmajor=True):
    '''
    Converts a csv file representing a numpy float64 (double) array into a binary file in row-major order
    :param inpath:
    :param outpath:
    :return:
    '''
    if rowmajor==False:
        indat=np.loadtxt(inpath,dtype=np.float64,delimiter=',')
    else:
        indat=np.loadtxt(inpath,dtype=np.float64,delimiter=',').transpose()
    outf=open(outpath,'wb')
    indat.tofile(outf)
    outf.close()
    del indat
    foin, fiin = os.path.split(inpath)
    foout, fiout = os.path.split(outpath)
    print ("Successfully converted %s to binary data, stored in file %s" % (fiin,fiout))

def convert_binary_to_csv(inpath,outpath,shape,rowmajor=True):
    '''
    Converts a csv file representing a numpy float64 (double) array into a binary file in row-major order
    :param inpath:
    :param outpath:
    :return:
    '''
    fin=open(inpath,'rb')
    arrnew=np.fromfile(fin,dtype=np.float64)
    fin.close()
    arr=np.reshape(arrnew,shape,'C')
    np.savetxt(outpath,arr,delimiter=',')


def test_numpy_inverse_timing(size):
    tic=dt.datetime.now()
    a=np.random.rand(size,size)
    b=np.linalg.inv(a)
    print ("time:\t%s" % str(dt.datetime.now()-tic))

# def make_kronecker_mvm_test_case(Aht, Awd,Bht,Bwd):
#     import cStringIO
#     a=np.random.rand(Aht,Awd)
#     b=np.random.rand(Bht,Bwd)
#     c=np.random.rand(Awd*Bwd,1)
#     f=cStringIO.StringIO()
#     f.write(np.array([a.shape[0]],dtype=np.int32).tobytes())
#     f.write(np.array([a.shape[1]],dtype=np.int32).tobytes())
#     f.write(np.array([b.shape[0]],dtype=np.int32).tobytes())
#     f.write(np.array([b.shape[1]],dtype=np.int32).tobytes())
#     f.write(a.tobytes())
#     f.write(b.tobytes())
#     f.write(c.tobytes())
#     g=open('test_binary_kronmvm.bin','wb')
#     g.write(f.getvalue())
#     g.close()
#     f.close()
#     a.tofile(open('test_a.bin','wb'))
#     b.tofile(open('test_b.bin','wb'))
#     c.tofile(open('test_c.bin','wb'))
#     np.savetxt('test_a.csv',a,delimiter=',')
#     np.savetxt('test_b.csv',b,delimiter=',')
#     np.savetxt('test_c.csv',c,delimiter=',')
#     out_shouldbe=np.dot(np.kron(a,b),c)
#     return a,b,c
#     print out_shouldbe

if __name__=='__main__':
    row_major=True

    # biom-data.csv:
    bdin=os.path.join(npdata,'biom-data.csv')
    bdout=os.path.join(npdata,'biomdata.bin')
    # convert_numpy_csv_to_binary(bdin,bdout,row_major)

    # dMatrix-full.csv:
    din=os.path.join(csvbackup,'dMatrix-full.csv')
    dout=os.path.join(npdata,'dmatrixfull.bin')
    # convert_numpy_csv_to_binary(din,dout,row_major)

    d=np.loadtxt(din,dtype=np.float64,delimiter=',')
    dL, dQ = np.linalg.eig(d)
    np.savetxt(os.path.join(csvbackup,'dMatrix-eigenQ.csv'),dQ,delimiter=',')
    np.savetxt(os.path.join(csvbackup,'dMatrix-eigenVals.csv'),dL,delimiter=',')

    f=open(os.path.join(npdata,'dmatrixEigenQ.bin'),'wb')
    dQ.tofile(f)
    f.close()
    f1=open(os.path.join(npdata,'dmatrixEigenVals.bin'),'wb')
    dL.tofile(f1)
    f1.close()
    # x-data.csv:
    xin=os.path.join(npdata,'x-data.csv')
    xout=os.path.join(npdata,'xdata.bin')
    # convert_numpy_csv_to_binary(xin,xout,row_major)

    # y-data.csv:
    yin=os.path.join(npdata,'y-data.csv')
    yout=os.path.join(npdata,'ydata.bin')
    # convert_numpy_csv_to_binary(yin,yout,row_major)
else:
    din=os.path.join(npdata,'dMatrix-full.csv')
    dinvin=os.path.join(npdata,'dMatrix-inverse-full.csv')
    betabin=os.path.join(npdata,'beta.bin')
    # d=np.loadtxt(din,dtype=np.float64,delimiter=',')
    # dinv=np.loadtxt(dinvin,dtype=np.float64,delimiter=',')
    # beta=np.fromfile(open(betabin,'rb'),dtype=np.float64)
