import pyomo.environ as pyo
from pyomo.opt import SolverFactory
import math
import random
import time

def time_reps(func, iteration_limit=100, time_limit=10):
    start = time.time()
    reps = 0
    for i in range(0, iteration_limit):
        func()
        reps += 1
        if time.time() - start > time_limit:
            break
    end = time.time()
    avg_ms = (end - start) / reps * 1000
    print("%s => %.3f ms" % (func.__name__, avg_ms))
    return

def perf_pyomo_sum():
    model = pyo.ConcreteModel()
    model.x = pyo.Var()
    model.obj = pyo.Objective(expr=sum(model.x**i for i in range(10000)))
    return

def perf_pyomo_prod():
    model = pyo.ConcreteModel()
    model.x = pyo.Var()
    model.obj = pyo.Objective(expr=math.prod(model.x**i for i in range(10000)))
    return

def perf_pyomo_many_constraints():
    model = pyo.ConcreteModel()
    model.X = pyo.RangeSet(0, 10000)
    model.x = pyo.Var(model.X)
    def constraint(model, i):
        return pyo.sin(model.x[i]) <= pyo.cos(i)
    model.c = pyo.Constraint(model.X, rule=constraint)
    return

def perf_pyomo_mle():
    model = pyo.ConcreteModel()
    n = 1000
    model.x = pyo.Var(initialize=0.0)
    model.y = pyo.Var(within=pyo.NonNegativeReals, initialize=1.0)
    data = [random.random() for _ in range(n)]
    model.obj = pyo.Objective(
        expr = n / 2 * pyo.log(1 / (2 * math.pi * model.y**2)) -
            sum((data[i] - model.x)**2 for i in range(n)) / (2 * model.y**2),
        sense = pyo.maximize,
    )
    opt = SolverFactory("ipopt")
    opt.solve(model, tee=False)
    return

def perf_pyomo_clnlbeam():
    N = 1000
    h = 1 / N
    alpha = 350
    model = pyo.ConcreteModel()
    model.S = pyo.RangeSet(1,N+1)
    model.S2 = pyo.RangeSet(1,N)
    model.t = pyo.Var(model.S, bounds=(-1.0, 1.0))
    model.x = pyo.Var(model.S, bounds=(-0.05, 0.05))
    model.u = pyo.Var(model.S)
    model.obj = pyo.Objective(
        expr = sum(
            0.5 * h * (model.u[i+1]**2 + model.u[i]**2) +
            0.5 * alpha * h * (pyo.cos(model.t[i+1]) + pyo.cos(model.t[i]))
            for i in model.S2
        )
    )
    def con_1(model, i):
        return model.x[i+1] - model.x[i] - 0.5 * h * (pyo.sin(model.t[i+1]) + pyo.sin(model.t[i])) == 0
    model.c1 = pyo.Constraint(model.S2, rule=con_1)
    def con_2(model, i):
        return model.t[i+1] - model.t[i] - 0.5 * h * model.u[i+1] - 0.5 * h * model.u[i] == 0
    model.c2 = pyo.Constraint(model.S2, rule=con_2)
    opt = SolverFactory("ipopt")
    opt.solve(model, tee=False)
    return

def perf_pyomo_rosenbrock():
    model = pyo.ConcreteModel()
    model.x = pyo.Var()
    model.y = pyo.Var()
    model.obj = pyo.Objective(
        expr = (1 - model.x)**2 + 100 * (model.y - model.x**2)**2
    )
    opt = SolverFactory("ipopt")
    opt.solve(model, tee=False)
    return

def perf_pyomo_jump_2788():
    N = 400
    k = N
    n = 12
    p = [random.randint(400, 700) for _ in range(k)]
    c1 = [[random.randint(100, 200) for _ in range(k)] for _ in range(n)]
    b = [random.randint(150, 250) for _ in range(k)]
    model = pyo.ConcreteModel()
    model.S = pyo.RangeSet(1, n)
    model.K = pyo.RangeSet(1, k)
    model.x = pyo.Var(model.S, bounds=(0, 1))
    model.var1 = pyo.Var(bounds=(0, 1))
    model.var2 = pyo.Var(bounds=(0, 1))
    model.var3 = pyo.Var(bounds=(0, 1))
    model.obj = pyo.Objective(
        expr=model.var1 - model.var2 + model.var3,
        sense=pyo.maximize,
    )
    model.expr = sum(model.x[i] * p[i-1] for i in model.S)
    def expr_c1(model, j):
        return sum(model.x[i] * c1[i-1][j-1] for i in model.S)
    def expr_c2(model, j):
        return sum(model.x[i] * 0.9 * c1[i-1][j-1] for i in model.S)
    model.con1 = pyo.Constraint(
        expr=model.expr==sum(b[j-1]/(1+model.var1)**j for j in model.K),
    )
    model.con2 = pyo.Constraint(
        expr=model.expr==sum(expr_c1(model, j)/(1+model.var2)**j for j in model.K),
    )
    model.con3 = pyo.Constraint(
        expr=model.expr==sum(expr_c2(model, j)/(1+model.var3)**j for j in model.K)
    )
    def con_4(model, j):
        return expr_c1(model, j) >= b[j-1]
    model.con4 = pyo.Constraint(model.K, rule=con_4)
    opt = SolverFactory("ipopt")
    opt.solve(model, tee=False)
    return

if __name__ == "__main__":
    for f in [
        perf_pyomo_sum, 
        perf_pyomo_prod, 
        perf_pyomo_many_constraints, 
        perf_pyomo_mle,
        perf_pyomo_clnlbeam,
        perf_pyomo_rosenbrock,
        perf_pyomo_jump_2788,
    ]:
        time_reps(f)
