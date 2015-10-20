import sympy as sym

chi = sym.Symbol('chi')
phi = sym.Symbol('phi')

a = sym.exp(1j*(chi+phi));
b = (sym.exp(1j*(2*chi)) + sym.exp(1j*(2*chi+2*phi)))/4;

c = sym.exp(1j*(chi-phi));
d = (sym.exp(1j*(2*chi)) + sym.exp(1j*(2*chi-2*phi)))/4;

tu = sym.Symbol('tu')
td = sym.Symbol('td')
 

M = sym.eye(2) - sym.Matrix([[b, -b],[-d,d]]);  
c = sym.Matrix([[a],[c]]);

t = M.inv() * c

t_total = sym.simplify( (t[0]+t[1])/2);
 

print(sym.latex(t_total))