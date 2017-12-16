
import numpy as np
import scipy

def Triangle_Algorithm(dataset, p, epsilon, **kwargs):
    matA = dataset
    [m, n] = matA.shape
    diffmat = np.subtract(matA, np.tile(p, (1, n)))

    eudis = np.sum(np.square(diffmat), axis=0)
    min_index = np.argmin(eudis)
    p_p = matA[:,int(min_index)]
    alpha = np.zeros((1, n), dtype=np.int)
    alpha[0][int(min_index)] = 1

    # Situation when input matrix has only one point
    distance = np.amax(eudis)
    if n<=1:
        if np.amin(eudis) != 0:
            return [0, p_p, alpha, np.amin(eudis)]
        else:
            return [1, p, alpha, np.amin(eudis)]

    # Iterative step
    inorout = 1
    iter = 0
    while float(np.sqrt(np.dot((p - p_p).conj().T, p - p_p))) > epsilon:
        iter = iter + 1
        found = 0
        distance = np.sqrt(np.dot((p-p_p).conj().T, p-p_p))
        rnd_index = np.random.choice(n, size=n, replace=False)
        for ii_index in rnd_index:
            iter = iter + 1
            if np.dot((p - matA[:,ii_index]).conj().T, p - matA[:,ii_index]) < np.dot((p_p - matA[:,ii_index]).conj().T, p_p - matA[:,ii_index]):
                beta = np.divide(float(np.dot((p_p-matA[:,ii_index]).conj().T, p_p-p)), float(np.dot((p_p-matA[:,ii_index]).conj().T, p_p-matA[:,ii_index])))
                alpha = np.dot(1 - beta, alpha)
                alpha[0][ii_index] = alpha[0][ii_index] + beta
                p_p = np.dot(1 - beta, p_p) + np.dot(beta, matA[:,ii_index])
                found = 1
                break

        if found == 0:
            inorout = 0
            break

    return [inorout, p_p, alpha, float(distance)]

# define input variables
# input_mat = np.matrix(np.arange(10).reshape((2,5)))
# point = np.matrix(np.arange(2).reshape((2,1)))
input_mat = np.matrix([[0,0,4,4],[0,4,0,4]])
point = np.matrix([[5],[5]])

# point = np.array([[0, 1]])
print 'Input matrix: '
print input_mat
print 'Point P: ' + str(np.matrix.flatten(point))

flag, p_prime, coef, distance = Triangle_Algorithm(dataset=input_mat, p=point, epsilon=0.1)

print '\nOutput: '

if flag == 0:
    print 'p is not in Convex Hull'
else:
    print 'p is in Convex Hull'

print 'p prime is at: ' + str(np.matrix.flatten(p_prime))

print 'distace of p from p prime is: ' + str(distance)
