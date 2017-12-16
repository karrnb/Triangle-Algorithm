
import numpy as np
import scipy

__author__ = 'karrnb'

def Triangle_Algorithm(data, p, epsilon, **kwargs):
    matA = data
    [m, n] = matA.shape
    diffmat = np.subtract(matA, np.tile(p, (1, n)))
    euclid_dist = np.sum(np.square(diffmat), axis=0)
    min_index = np.argmin(euclid_dist)
    p_prime = matA[:,int(min_index)]
    alpha = np.zeros((1, n), dtype=np.int)
    alpha[0][int(min_index)] = 1

    # Situation when input matrix has only one point
    distance = np.amin(euclid_dist)
    if n<=1:
        if np.amin(euclid_dist) != 0:
            return [0, p_prime, alpha, np.amin(euclid_dist)]
        else:
            return [1, p, alpha, np.amin(euclid_dist)]

    # Iterative step
    flag = 1
    iter = 0
    while float(np.sqrt(np.dot((p - p_prime).conj().T, p - p_prime))) > epsilon:
        iter = iter + 1
        found = 0
        distance = np.sqrt(np.dot((p - p_prime).conj().T, p - p_prime))
        random_index = np.random.choice(n, size=n, replace=False)
        for index in random_index:
            iter = iter + 1
            if np.dot((p - matA[:,index]).conj().T, p - matA[:,index]) < np.dot((p_prime - matA[:,index]).conj().T, p_prime - matA[:,index]):
                beta = np.divide(float(np.dot((p_prime - matA[:,index]).conj().T, p_prime - p)), float(np.dot((p_prime - matA[:,index]).conj().T, p_prime - matA[:,index])))
                alpha = np.dot(1 - beta, alpha)
                alpha[0][index] = alpha[0][index] + beta
                p_prime = np.dot(1 - beta, p_prime) + np.dot(beta, matA[:,index])
                found = 1
                break

        if found == 0:
            flag = 0
            break

    return [flag, p_prime, alpha, float(distance)]

# define input variables
# input_mat = np.matrix(np.arange(10).reshape((2,5)))
# point = np.matrix(np.arange(2).reshape((2,1)))
input_mat = np.matrix([[0,0,4,4],[0,4,0,4]])
point = np.matrix([[5],[5]])

# point = np.array([[0, 1]])
print 'Input matrix: '
print input_mat
print 'Point P: ' + str(np.matrix.flatten(point))

flag, p_primerime, coef, distance = Triangle_Algorithm(data=input_mat, p=point, epsilon=0.1)

print '\nOutput: '

if flag == 0:
    print 'p is not in Convex Hull'
else:
    print 'p is in Convex Hull'

print 'p prime is at: ' + str(np.matrix.flatten(p_primerime))

print 'distace of p from p prime is: ' + str(distance)
