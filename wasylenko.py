# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

def vecAvg(vec):
    '''Adds up all the values in a numpy array, and divide by the count of values.
    Assumes that vec is a numpy ndarray or numpy matrix.'''
    return np.sum(vec) / vec.size

def doint(initial, dxdt, tmin=0, tmax=800, nstep=1e5):
    '''Integrate a given differential system forward from given initial conditions.
    Returns X and T:
     in X, the first index is across times, and the second index is across state variables.
     T has only one index; across times.
    It's possible that Scipy is using adaptive time stepping, in which case the (linearly
     spaced) values of T will be inaccurate. Use with caution.'''
    def dxdtextra(X, T=0):
        return dxdt(X, T=T)
    T = np.linspace(tmin, tmax, nstep)
    X0 = np.array(initial)
    X, infodict = integrate.odeint(dxdtextra, X0, T, full_output=True)
    return (X, T)

def main(N=60, omega=6, s=.5, plot=True, transient=True, tmax=.125, nstep=1000, verbose=True, X0=None):
    
    # parameters and auxiliary functions:
    
       # These are intended as vector-valued functions--they take a vector,
       #   and return a vector (i.e., a numpy ndarray or matrix).
       #   However, they would also work just fine with scalars. Python is duck-typed.
    def minf(v):
        return 1 / (1 + np.exp(-(v - thm) / sigm))
    
    def hinf(v):
        return 1 / (1 + np.exp(-(v - thh) / sigh))
    
    def sinf(v):
        return 1 / (1 + np.exp(-(v - ths) / sigs))
        
    def tauinf(v):
        return tau1 + (tau2 - tau1) / (1 + np.exp(-(v - thtau) / sigtau))
    
    # fixed parameters
    ths = -20
    sigs = 2
    thh = -79
    sigh = -5
    thm = -65
    sigm = 7.8
    thtau = -65
    sigtau = 4
    tau1 = 1
    tau2 = 80
    
    gLTC = .01
    VLTC = -75
    VsynTC = 0
    gLRE = .2
    VLRE = -80
    epRE = 2
    VsynRE = -80
    gCa = 1
    Cm = 1
    
    # bifurcation parameters
    epTC = 1 + 2 * s
    gTC = 0.03 + 0.07 * s
    gRE = 0.1 + 0.2 * s
    VCa = 120 - 30 * s
    
    def theSum(sinfVTC):
        '''Used for computing the sum of inputs to each of the
        RE neurons. Takes as argument a precomputed vector of s_inf(v_i)
        values for each of the TC neurons' voltages, v_i.'''
        shape = sinfVTC.shape
        sinfVTC = list(sinfVTC)
        N = len(sinfVTC)
        toreturn = list(np.zeros(shape))
        for i in range(N):
#             print "---"
            for k in range(-omega, omega+1):
                index = i + k
                if index >= N: index -= N
#                 print "i:", i, "k:", k,  "index:", index
                toreturn[i] += sinfVTC[index]
        return np.array(toreturn).reshape(shape)
    
#     wuts = [[], [], []]
    def diff(X, T=0):
        '''The vector system-of-equations.
        For the state vector X, of shape (4*N,)
            the first  N elements are the VTC voltages,
            the second N elements are the hTC coeffecients
            the third  N elements are the VRE voltages, and
            the fourth N elements are the hRE coeffecients.
        Returns an array of shape (4*N,) of the rates-of-change of these state variables.
        '''
        # split the state vector into meaningful pieces:
        VTC = X[0*N:1*N]
        hTC = X[1*N:2*N]
        VRE = X[2*N:3*N]
        hRE = X[3*N:4*N]
        sinfVTC = sinf(VTC)  # precompute s_inf(v) values for use in theSum()
        
        # the four differential equations:
        dVTCdt = (-gLTC * (VTC - VLTC)
                  -gCa * (minf(VTC)) ** 3 * hTC * (VTC - VCa)
                  -gTC * sinf(VRE) * (VTC - VsynRE)
                 ) / Cm
        dhTCdt = epTC * (hinf(VTC) - hTC) / tauinf(VTC)
        dVREdt = (-gLRE * (VRE - VLRE)
                  -gCa * (minf(VRE)) ** 3 * hRE * (VRE - VCa)
                  -gRE / (2*omega+1) * theSum(sinfVTC) * (VRE - VsynTC)
                 ) / Cm
        dhREdt = epRE * (hinf(VRE) - hRE) / tauinf(VRE)
        # show progress
        if verbose:
            if not int(T) % 100:
                print T
        
        # Wasylenko 2010, Sec. 3.2:
        #  First, we started with a perturbation in the vTC from the (spatially
        #  uniform) rest state. This generated both left- and right-traveling 
        #  waves, but by transiently changing the boundary conditions it was 
        #  possible to eliminate the left-travelingwave.
        # Here is that transient modification. It's possible that I could do it better:
        if transient:
            if T<600:
                dVTCdt[0:omega] = 0
                dhTCdt[0:omega] = 0
                dVREdt[0:omega] = 0
                dhREdt[0:omega] = 0
            
#         if 600 < T < 610:
#             pp(X)
            
        return np.array(np.vstack((dVTCdt, dhTCdt, dVREdt, dhREdt))).reshape((4*N,))

    # Initial Conditions (near resting point,
    #  which was determined by poorly choosing
    #  ICs and then measuring their mean steady-state values.)
    if transient:
        X0 = np.ones((N * 4,))
        X0[0*N:1*N] = -49  # VTC_0
        X0[1*N:2*N] = 0    # hTC_0
        X0[2*N:3*N] = -78  # VRE_0
        X0[3*N:4*N] = 0.5  # hRE_0
        X0[omega+3] = 0    # a small perturbation in VTC_0
    
    # Uncomment this for krazy fun. Surprisingly, the system still settles
    #  into a traveling-wave behavior, if the transient BC-modification is maintained long enough:
#     X0 = np.random.random((N*4,)) * 120 - 60 

    X, T = doint(X0, diff, tmin=0, tmax=tmax, nstep=nstep)
    
    if plot:
        plotMultiple(X, T, s=s)
        
    # This is what I used to determine the "resting" conditions
    #  (except, without the transient BC-modification):
    if verbose:
        print "average final VTC is", vecAvg(X[-1, 0*N:1*N])
        print "average final hTC is", vecAvg(X[-1, 1*N:2*N])
        print "average final VRE is", vecAvg(X[-1, 2*N:3*N])
        print "average final hRE is", vecAvg(X[-1, 3*N:4*N])
    return X, T

def plotMultiple(X, T, s=0.5, show=False):
    N = X.shape[1] / 4
    N = int(N)
    fig = plt.figure(figsize=(11*1.5, 8.5*1.5))
    axes = [fig.add_subplot(1, 4, i+1) for i in range(4)]
    imgplots = range(4)
    cmaps = ['hot', 'jet', 'cool', 'bone']
    titles = [r"$V^{TC}$", r"$h^{TC}$", r"$V^{RE}$", r"$h^{RE}$"]
    for i in range(4):
        imgplots[i] = axes[i].imshow(X[:, i*N:(i+1)*N], aspect="auto", origin="lower")
        axes[i].set_title(titles[i])
        imgplots[i].set_interpolation('nearest')
        imgplots[i].set_cmap(cmaps[i])
        fig.colorbar(imgplots[i], orientation ='horizontal')
    fig.suptitle(r"$s=%.2f$" % s)
    fig.savefig("wasylenko_waves-s%.2f.pdf" % s)

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(1, 1, 1)
    for i in range(N):
        ax2.plot(X[-50:-1, i], X[-50:-1, 2*N+i], 'k')
    ax2.set_xlabel(r"$V^{TC}$")
    ax2.set_ylabel(r"$V^{RE}$")
    fig2.suptitle(r"$s=%.4f$" % s)
    fig2.savefig("phasespace-s%.4f.pdf" % s)

    fig3 = plt.figure()
    ax3 = fig3.gca()
    im = ax3.imshow(X[:, 0:N], aspect="auto", origin="lower")
    im.set_interpolation('nearest')
    im.set_extent([0, N, min(T), max(T)])
    ax3.set_xlabel(r'$V^{TC}$ neuron index')
    ax3.set_ylabel('time')
    fig3.colorbar(im)

    if show:
        plt.show()

def sRange():
    '''Run the integration for several different values of s'''
    fig4 = plt.figure()
    ax4 = fig4.add_subplot(1, 1, 1)
    s = list(np.arange(0, 0.5, 0.1))
    s.extend(list(np.arange(0.5, 1, 0.02)))
    avg = range(len(s))
    for i, s_ in enumerate(s):
        X, T = main(N=60, s=s_, omega=6, tmax=3000, nstep=2000)
        avg[i] = np.sum(X[:, 0:N]) / X[:, 0:N].size
    ax4.plot(s, avg)
    ax4.set_xlabel(r"arc length parameter, $s$")
    ax4.set_ylabel(r"average $V^{TC}$")
    fig4.suptitle("$V^{TC}$ averaged across %d neurons for different values of $s$" % N)
    fig4.savefig("averge_VTC_values_over_s.pdf")

def poi(s, tmax=2000, N=60, omega=6):
    # trying to do fig 4
    X, T = main(tmax=tmax, s=s, N=N, omega=omega)
    #N = X.shape[1]
    N = 20
    Tmax = X.shape[0]
    fig = plt.figure()
    VTC = X[0:N,-1]
    ax = fig.gca()
    #for j, color in zip([600, 700, 800], ['red', 'blue', 'green']):
    color='green'
    #for j, color in zip([800], ['green']):
    cutIndices = []
    for p0 in range(-10, 50):
        for j in range(585, 1000):
            for i in range(N):
                if  p0 < X[j,i] < p0+1:  # take cuts near-ish to peak
                    cutIndices.append((j, i))
    vp = []
    vm = []
    for cutIndex in cutIndices:
        left = list(cutIndex)
        left[1] -= 1
        right = list(cutIndex)
        right[1] += 1
        
        if left[1] < 0:
            left[1] += N
        if right[1] >= N:
            right[1] -= N
        print cutIndex, left, right
        try:
            vp.append(X[left])
            vm.append(X[right])
            #ax.scatter([i], [j])
        except:
            pass
    ax.scatter(vp, vm)
    #ax.legend(loc="best")
    ax.set_xlabel(r"$V^{TC}_{j+1}$")
    ax.set_ylabel(r"$V^{TC}_{j-1}$")
    #ax.plot(X[800, 12:34])
    plt.show()
    
def pMap(X0, tau=None, d=1, s=0.8, omega=6, N=60):
    if tau==None:
        tau = findTau(X0, d=d, N=N, omega=omega, s=s)
    X0 = shift(X0, d=d)
    Xtau, T = main(N=N, omega=omega, s=s, plot=False, tmax=tau, verbose=False, X0=X0)
    return Xtau[-1,:].reshape(X0.shape)

def shift(X, d=1):
    if d==0:
        return X
    #if d<0:
        #print "d<0 not implemented."
    originalShape = X.shape
    X = X.reshape((X.size, ))
    if abs(d)>1:
        #There's probably a better way to do this than with recursion,
        # but I'm feeling lazy right now.
        d = (abs(d) - 1) * d / abs(d)
        X = shift(X, d=d)
    if d<0:
        return np.hstack((X[-1], X[:-1])).reshape(originalShape)
    else:
        return np.hstack((X[1:], X[0])).reshape(originalShape)
    
def shiftSolnVec(X, d=1):
    originalShape = X.shape
    X = X.reshape((X.size, ))
    N = X.size / 4
    Vl = [X[i*N:(i+1)*N] for i in range(4)]
    Vl = [shift(V, d=d) for V in Vl]
    return np.hstack(Vl).reshape(originalShape)
    
def shiftSolnArray(A, d=1):
    for t in range(A.shape[0]):
        A[t,:] = shiftSolnVec(A[t,:], d=d)
    return A

def findTau(X0, d=1, pickPoint=0.5, testValue=None, searchTime=650, N=60, omega=6, s=0.8):
    originalShape = X0.shape
    X0 = X0.reshape((X0.size, ))
    testIndex = int(X0.size * pickPoint)
    #testIndex = X0.argmax()
    if testValue==None:
        testValue = X0[testIndex]
    if d==0:
        return 0
    X0 = shiftSolnVec(X0, d=d)
    print "start at", X0[testIndex], ", search for", testValue
    def integrator(X0, tmax=searchTime):
        return main(N=N, omega=omega, s=s, plot=False, transient=False, tmax=tmax, verbose=False, X0=X0)
    Xf, T = integrator(X0)
    closestDiff = min(list(Xf[:,testIndex] - testValue))
    count = 0
    #plotMultiple(Xf, T, show=True)
    while abs(closestDiff) > 10:
        count += 1
        print count
        Xf, T = integrator(Xf[-1, :], tmax=searchTime)
        closestDiff = min(list(Xf[:,testIndex] - testValue))
        closestIndex = list(Xf[:,testIndex] - testValue).index(closestDiff)
        print "closest is", Xf[closestIndex, testIndex], "with diff", closestDiff
        if count > 1:
            print "probably looking the wrong way. try h-flipping data"
        if count > 4:
            print "this is taking too long. breaking."
            break
    closestIndex = list(Xf[:,testIndex] - testValue).index(closestDiff)
    print closestDiff, closestIndex, T.min(), "<", T[closestIndex], "<", T.max()
    return searchTime * count + T[closestIndex]
    
def testPMap():
    from wasyData import loadData
    data = loadData('wasys8.dat')
    X = data[:, 1:]
    T = data[:, 0]
    #plotMultiple(X, T, s=0.5, show=True)
    d = 1
    #plotMultiple(shiftSolnArray(X, d=d), T, show=True)
    X0 = data[60][1:]
    #whatThis(X0=X0)
    print findTau(X0, d=d)

def whatThis(X0=None):
    if X0==None:
        from wasyData import loadData
        X0 = loadData('wasys8.dat')[60][1:]
    N = X0.size / 4
    f = plt.figure()
    a = f.gca()
    a.plot(X0[0:N])
    #for tmax in [200, 400, 600]:
        #X, T = main(tmax=tmax, X0=X0, plot=False)
        #a.plot(X[-1,:][0:50])
    plt.show()
    

if __name__=="__main__":
    #main(N=60, omega=6, tmax=2000, s=0.6, nstep=1000)
    # uncomment either the sRange line to run for many values of s:
#     sRange()
    # In that case, you might want to comment out the plt.show line.
    #poi(0.715)
    #plt.show()
    testPMap()
    #whatThis()

