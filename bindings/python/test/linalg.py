import moopy.linalg as linalg

def dot():
    return linalg.dot([1, 2], [3, 4]) == 11

def square():
    x = [1.0, 2.0]
    x = linalg.square(x)
    return (x[0] == 1.0) and (x[1] == 4.0)

def matrix_vector():
    A = [1.0, 2.0, 3.0, 4.0]
    x = [5.0, 6.0]
    b = linalg.matrix_vector('R', A, x)
    return (b[0] == 17.0) and (b[1] == 39.0)    

def diagmat_vec():
    D = [1.0, 2.0]
    x = [4.0, 5.0]
    b = linalg.diagmat_vec(D, False, x)
    return (b[0] == 4.0) and (b[1] == 10.0)

def dsaxpy():
    x = [1.0, 2.0]
    beta = 2.0
    y = [4.0, 5.0]
    D = [3.0, -3.0]
    y_out = linalg.dsaxpy(x, y, D, beta, False)
    return (y_out[0] == 27.0) and (y_out[1] == -36.0)

def dgmv():
    x = [5.0, 6.0]
    y = [5.0, 6.0]
    D = [7.0, 8.0]
    beta = -2.0
    b = linalg.dgmv(x, y, D, beta, False)
    return (b[0] == 25.0) and (b[1] == 36.0)

print(dot())
print(square())
print(matrix_vector())
print(diagmat_vec())
print(dsaxpy())
print(dgmv())
