import sympy as s

rho, xi, eta, chi = s.symbols('rho, xi, eta, chi')
mp = 4 * rho * (1 + etap) / ((1 + eta + rho)**2 + (xi - chi)**2)
K = s.Function('K')
phi = K(mp) * s.sqrt((1 + eta) * mp / rho)
expr = phi.diff(chi, 2)
sol = expr.subs([(chi, 0), (eta, 0)])
print(s.latex(sol))