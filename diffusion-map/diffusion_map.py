import scipy.linalg as la
import numpy as np
import matplotlib.pyplot as plt

def ring_matrix(n,t):
    first_row = np.zeros(n)
    first_row[1]=1
    first_row[n-1]=1
    W = np.matrix(la.circulant(first_row))
    D = np.matrix(np.diag(np.ones(n)*2))
    D_inverse = la.inv(D)
    M = np.dot(D_inverse, W)    
    D_sqrt = np.sqrt(D)
    D_sqrt_inverse = la.inv(D_sqrt)
    S = np.dot(np.dot(D_sqrt_inverse, W), D_sqrt_inverse)
    dftmtx = np.matrix(np.fft.fft(np.eye(n)))/np.sqrt(n)
    Lambda = np.dot(np.dot(dftmtx.getH(), S), dftmtx)
    if n % 2 == 0:
        index_1 = int((n-2)/2)
        index_2 = int((n+2)/2)
    else:
        index_1 = int((n-1)/2)
        index_2 = int((n+1)/2)
    lambda_1 = Lambda[index_1,index_1].real
    lambda_2 = Lambda[index_2,index_2].real
    v_1 = dftmtx[index_1]
    v_2 = dftmtx[index_2]
    u_1 = ((v_1+v_2)/2).real
    u_2 = (1j*(v_1-v_2)/2).real
    p_1 = u_1/np.linalg.norm(u_1)
    p_2_temp = u_2 - float(u_2*p_1.transpose())*p_1
    p_2 = p_2_temp/np.linalg.norm(p_2_temp)
    phi_1 = p_1/2
    phi_2 = p_2/2
    row_dim1 = lambda_1 ** t * phi_1
    row_dim2 = lambda_2 ** t * phi_2
    x,y = row_dim1.tolist()[0],row_dim2.tolist()[0]
    plt.plot(x,y,"o")
    plt.axis('equal')
    plt.show()

ring_matrix(7,10)
ring_matrix(8,10)