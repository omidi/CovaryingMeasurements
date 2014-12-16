

import numpy as np


def inv(X):
    return np.linalg.inv(X)


def T(X):
    return np.transpose(X)


def det(X):
    return np.linalg.det(X)


def prod(X,Y,Z):
    return np.dot(np.dot(X, Y), Z)


def dependencyScore(X, Y, C, S):
    n, c = X.shape
    C_inv = inv(C)
    S_inv = inv(S)
    N = inv(2. * C + np.dot(np.dot(C, S_inv), C))
    M = inv(C + np.dot(np.dot(C, S_inv), C))
    XY = X + Y
    ###
    beta = 0. 
    for xy in XY.T:
        beta += .5 * prod( xy ,N ,xy.T )
    for x in X.T:
        beta -= .5 * prod( x, M, x.T )
    for y in Y.T:
        beta -= .5 * prod( y, M, y.T )
    temp = .5*np.log( det(S) )
    temp += np.log( det( inv(2.0*C_inv + S_inv)) )
    temp -= 2.*np.log( det(inv(C_inv + S_inv)) )
    beta += 0.5 * c * temp
    return beta


def flatLikelihoodScore(X, Y, C, S):
    n, c = X.shape
    C_inv = inv(C)
    S_inv = inv(S)
    N = inv(2. * C + np.dot(np.dot(C, S_inv), C))
    K = inv( 2.* c *C + np.dot(np.dot(C, S_inv), C))
    XY = X + Y
    ###
    beta = 0.
    tmp = np.zeros(n)
    for xy in XY.T:
        beta += .5 * prod( xy ,N ,xy.T )
        tmp += xy
    beta -= .5 * prod( tmp ,K ,tmp.T )
    beta += .5 * (1. - c) * np.log( det(S) )
    beta += .5 * c * np.log( det(C) )
    beta -= .5 * np.log( det( inv(2.0*c*C_inv + S_inv)) )    
    beta += .5 * c * np.log( det( inv(2.0*C_inv + S_inv)) )
    return beta



def constantFunction(X, Y, C, S):
    c, n = X.shape
    C_inv = inv(C)
    S_inv = inv(S)
    XY = X+Y
    # M = inv(C + np.dot(np.dot(C, S_inv), C))
    N = inv(2.*c*C + np.dot(np.dot(C, S_inv), C))
    likelihood = 0.
    likelihood -= (c+1)*n*np.log(2. * np.pi)
    likelihood -= .5*( np.log(det(C)) - np.log(det(S)) - np.log(det((2*c*C_inv + S_inv))) )
    for x in X.T:
        likelihood -= .5 * prod( x, C_inv, x.T )
    for y in Y.T:
        likelihood -= .5 * prod( y, C_inv, y.T )
    for xy in XY.T:
        likelihood += .5 * prod(xy, N, xy.T)
    likelihood += n*np.log( 2. * np.pi )
    return likelihood



def marginalLikelihood(X, C, S):
    c, n = X.shape
    C_inv = inv(C)
    S_inv = inv(S)
    M = inv(C + np.dot(np.dot(C, S_inv), C))
    likelihood = 0.
    for x in X.T:
        likelihood += .5 * prod( x, M, x.T )
        likelihood -= .5 * prod( x, C_inv, x.T )
        # print .5 * prod( x, M, x.T ) - .5 * prod( x, C_inv, x.T )
    likelihood -= .5*c*n*np.log(2. * np.pi)
    likelihood -= .5*c*(np.log(det(C)) + np.log(det(S)) - np.log(det(C_inv + S_inv)))
    return likelihood


def jointLikelihood(X, Y, C, S):
    c, n = X.shape
    C_inv = inv(C)
    S_inv = inv(S)
    XY = X+Y
    M = inv(C + np.dot(np.dot(C, S_inv), C))
    N = inv(2.*C + np.dot(np.dot(C, S_inv), C))
    likelihood = 0.
    for xy in XY.T:
        likelihood += .5 * prod(xy, N, xy.T)
    for x in X.T:
        likelihood -= .5 * prod( x, M, x.T )
    for y in Y.T:
        likelihood -= .5 * prod( y, M, y.T )
    likelihood -= c*n*np.log(2. * np.pi)
    likelihood -= .5*c*(np.log(det(C)) + np.log(det(S)) - np.log(det((2*C_inv + S_inv))))
    return likelihood


def convertToSlope(X):
    return (X[:, 1:] - X[:, :-1])
    

def main():
    c = 10
    n = 5
    scale= 1.
    # sigma = scale* np.array([.02, .01, .02, .012, .078, .0068, .0221, .00865, .092, .102])
    # sigma = scale* np.array([.000012, .000005, .00002, 0.00002, .0000132, .0000028, 0.00000221, 0.00001865, .0000092, .000012])
    # sigma = scale* np.array([1.150012, .900005, .82, .51])
    sigma = 10*np.array([0.05858489, 0.04778119 , 0.08578143, 0.03702142, 0.044235])
    # sigma = np.array([0.00005858489, 0.00004778119 , 0.00008578143, 0.00003702142])
    # sigma = np.array([5.10858489, 5.14778119 , 5.12578143, 5.23702142])
    for k in xrange(100):
        tmpX = []
        tmpY = []
        for v in sigma:
            tmpX.append(np.random.normal(0., 1., c))
            tmpY.append(np.random.normal(0., 1., c))
            # tmpX.append(np.random.rand(c) * 2.)
            # tmpY.append(np.random.rand(c) * 2.)
        # tmpX = np.repeat(np.random.normal(0,1,n), c)
        # tmpY = np.repeat(np.random.normal(0,1,n), c)
        X = np.array(tmpX).reshape(n, c)
        Y = np.array(tmpY).reshape(n, c)        
        Z = np.copy(X)
        Z = scale * Z
        F = convertToSlope(X)
        G = convertToSlope(Z)
        K = convertToSlope(Y)
        # adding the noise
        for v, i in zip(sigma, xrange(len(sigma))):
            X[i] += np.random.normal(0., v)
            Y[i] += np.random.normal(0., v)
            Z[i] += np.random.normal(0., v)        

        # X = (X- np.mean(X))/np.std(X)
        # Y = (Y- np.mean(Y))/np.std(Y)
        # Z = (Z- np.mean(Z))/np.std(Z)        
        C = sigma*np.identity(n)
        C_inv = inv(C)
        # S = sigma*np.identity(n)
        S = np.identity(n)
        S_inv = inv(S)
        X_likelihood = marginalLikelihood(X, C, S)
        Y_likelihood = marginalLikelihood(Y, C, S)
        Z_likelihood = marginalLikelihood(Z, C, S)
        XY_likelihood = jointLikelihood(X, Y, C, S)
        XX_likelihood = jointLikelihood(X, X, C, S)
        ZZ_likelihood = jointLikelihood(Z, Z, C, S)
        YZ_likelihood = jointLikelihood(Y, Z, C, S)
        XZ_likelihood = jointLikelihood(X, Z, C, S)
        YY_likelihood = jointLikelihood(Y, Y, C, S)
        
        # print XY_likelihood, XZ_likelihood, YZ_likelihood, X_likelihood, Y_likelihood, Z_likelihood, XX_likelihood, YY_likelihood
        # print XY_likelihood - X_likelihood - Y_likelihood
        # print XZ_likelihood - X_likelihood - Z_likelihood
        # print XX_likelihood - 2*X_likelihood
        # print YZ_likelihood - Y_likelihood - Z_likelihood
        # print dependencyScore(X, Y, C, S)
        # print dependencyScore(X, Y, C, S), '\t', dependencyScore(X, 1.*Z, C, S), '\t',dependencyScore(X, 2.*Z, C, S), '\t',dependencyScore(X, 3.*Z, C, S), '\t',dependencyScore(3*X, 0.9*Z, C, S), '\t',dependencyScore(3*X, Z, C, S), '\t', dependencyScore(X, X, C, S)
        print dependencyScore(X, Y, C, S), '\t', dependencyScore(X, Z, C, S), '\t', dependencyScore(X, X, C, S), '\t',  \
          flatLikelihoodScore(X, Y, C, S), '\t', dependencyScore(F, K, C, S), '\t', dependencyScore(F, G, C, S)
        # print XY_likelihood - X_likelihood - Y_likelihood, '\t', \
        #     XZ_likelihood - X_likelihood - Z_likelihood
    exit()
    N = inv(2. * C + np.dot(np.dot(C, S_inv), C))
    M = inv(C + np.dot(np.dot(C, S_inv), C))
    XY = X + Y
    ###
    beta = 0.
    for xy in XY.T:
        beta += .5 * prod( xy ,N ,xy.T )
    for x in X.T:
        beta -= .5 * prod( x, M, x.T )
    for y in Y.T: 
        beta -= .5 * prod( y, M, y.T )
    
    temp = np.log( det(C) ) + np.log( det(S) )
    temp -= np.log( det(inv(2.0*C_inv + S_inv)) )
    temp += 2*np.log( det(inv(C_inv + S_inv)) )
    beta += float(c)/2 * temp
    beta += float(n+c) / 2. * np.log(2* np.pi)
    with open("data", "w") as outf:
        for x, y in zip(X.T,Y.T):
            outf.write("%f\t%f\n" % (x[0],y[0]))

    X = np.array(tmpX).reshape(n, c)
    Y = X 
    C = sigma*np.identity(n)
    C_inv = inv(C)
    S = sigma*np.identity(n)
    S_inv = inv(S) 
    N = inv(2. * C + np.dot(np.dot(C, S_inv), C))
    M = inv(C + np.dot(np.dot(C, S_inv), C))
    XY = X + Y
    ###
    beta = 0. 
    for xy in XY.T:
        beta += .5 * prod( xy ,N ,xy.T )
    for x in X.T:
        beta -= .5 * prod( x, M, x.T )
    for y in Y.T: 
        beta -= .5 * prod( y, M, y.T )
    
    temp = np.log( det(C) ) + np.log( det(S) )
    temp -= np.log( det(inv(2.0*C_inv + S_inv)) )
    temp += 2*np.log( det(inv(C_inv + S_inv)) )
    beta += float(c)/2 * temp
    beta += float(n+c) / 2. * np.log(2* np.pi)
    print np.sum(beta) 
    
    
if __name__ == '__main__':
    main()
