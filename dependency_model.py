

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
    c, n = X.shape
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
    temp = np.log( det(C) ) + np.log( det(S) )
    temp += np.log( det( inv(2.0*C_inv + S_inv)) )
    temp -= np.log( det(inv(C_inv + S_inv)) )
    beta += float(c)/2 * temp
    # beta -= .5*float(n+c) * np.log(2* np.pi)
    return beta 
    

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


def main():
    c = 7
    n = 4
    scale= 1.
    # sigma = scale* np.array([.02, .01, .02, .012, .078, .0068, .0221, .00865, .092, .102])
    # sigma = scale* np.array([.000012, .000005, .00002, 0.00002, .0000132, .0000028, 0.00000221, 0.00001865, .0000092, .000012])
    # sigma = scale* np.array([1.150012, .900005, .82, .51])
    sigma = scale* np.array([.0901000012, .091000005, .20600002, .30800001])
    for k in xrange(100):
        tmpX = []
        tmpY = []
        for v in sigma:
            tmpX.append(np.random.normal(0., 1., c))
            tmpY.append(np.random.normal(0., 1., c))

        X = np.array(tmpX).reshape(n, c)
        Y = np.array(tmpY).reshape(n, c)
        Z = np.copy(X)
        # adding the noise
        for v, i in zip(sigma, xrange(len(sigma))):
            X[i] += np.random.normal(0., v)
            Y[i] += np.random.normal(0., v)
            Z[i] += np.random.normal(0., v)

        C = sigma*np.identity(n)
        C_inv = inv(C)
        S = sigma*np.identity(n)
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
        print dependencyScore(X, Y, C, S), '\t', dependencyScore(X, Z, C, S), '\t', dependencyScore(X, X, C, S)
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
    print np.sum(beta) 

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
