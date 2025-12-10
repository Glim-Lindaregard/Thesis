from tsl_optimizer.parameters import *  
import quadprog 
import cvxopt



class Solver : 
    def __init__(self): 
        self.P = P 
        self.q = q 
        self.G = G 
        self.h = h 
        self.A = A

    def __quadprog_solve_qp(P, q, G, h, A=None, b=None):
        qp_G = .5 * (P + P.T)   # make sure P is symmetric
        qp_a = -q
        if A is not None:
            qp_C = -np.vstack([A, G]).T
            qp_b = -np.hstack([b, h])
            meq = A.shape[0]
        else:  # no equality constraint
            qp_C = -G.T
            qp_b = -h
            meq = 0
        return quadprog.solve_qp(qp_G, qp_a, qp_C, qp_b, meq)[0]

    def __cvxopt_solve_qp(self, P, q, G, h, A=None, b=None):
        P = .5 * (P + P.T)  # make sure P is symmetric
        args = [cvxopt.matrix(P), cvxopt.matrix(q)]
        args.extend([cvxopt.matrix(G), cvxopt.matrix(h)])
        if A is not None:
            args.extend([cvxopt.matrix(A), cvxopt.matrix(b)])
        sol = cvxopt.solvers.qp(*args)
        if 'optimal' not in sol['status']:
            return None
        return np.array(sol['x']).reshape((P.shape[1],))  

    def solve(self, U): 
        return self.__cvxopt_solve_qp(self.P, self.q, self.G, self.h, self.A, U)
