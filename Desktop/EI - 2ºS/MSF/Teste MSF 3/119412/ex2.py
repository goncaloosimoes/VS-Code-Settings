import numpy as np
import matplotlib.pyplot as plt

# Constantes

MASSA = 2
K = 0.5
B = 0.2
OMEGA_F = 1.0
OMEGA_0 = np.sqrt(K / MASSA)
F0 = 5

# Resolução

def maxminv(xm1,xm2,xm3,ym1,ym2,ym3):  
    # Máximo ou mínimo usando o polinómio de Lagrange
    # Dados (input): (x0,y0), (x1,y1) e (x2,y2) 
    # Resultados (output): xm, ymax 
    xab=xm1-xm2
    xac=xm1-xm3
    xbc=xm2-xm3

    a=ym1/(xab*xac)
    b=-ym2/(xab*xbc)
    c=ym3/(xac*xbc)

    xmla=(b+c)*xm1+(a+c)*xm2+(a+b)*xm3
    xm=0.5*xmla/(a+b+c)

    xta=xm-xm1
    xtb=xm-xm2
    xtc=xm-xm3

    ymax=a*xtb*xtc+b*xta*xtc+c*xta*xtb
    return xm, ymax


def forca_externa(t, omega_f = OMEGA_F, f0 = F0):
    return f0 * np.cos(omega_f * t)

def forca(x, k = K, b = B):
    return - k * x

def forca_amortecimento(v, b = B):
    return - b * v

def energia_potencial(x, k = K):
    return 0.5 * k * (x**2)

def f(t, x, v, m = MASSA):
    return (forca(x) + forca_amortecimento(v) + forca_externa(t)) / m

def euler(x0, v0, tf, dt, f):
    n = int(np.ceil(tf/dt))
    a = np.zeros(n + 1)
    v = np.zeros(n + 1)
    x = np.zeros(n + 1)
    t = np.zeros(n + 1)

    x[0] = x0
    v[0] = v0

    for i in range(n):
        a[i + 1] = f(t[i], x[i], v[i])

        v[i + 1] = v[i] + a[i] * dt
        x[i + 1] = x[i] + v[i] * dt
        t[i + 1] = t[i] + dt
    
    return a, v, x, t

x_eq = 0
x0 = 0.2
v0 = 0
tf = 200
dt = 0.01

# Alínea A
a, v, x, t = euler(x0, v0, tf, dt, f)

pontos_criticos = []
for i, vel in enumerate(v):
    if i < v.size - 1 and v[i] * v[i+1] < 0:
        pontos_criticos.append(i)

pontos_criticos = np.array(pontos_criticos)

t_crit = []
x_crit = []

for i in pontos_criticos:
    _t, _x = maxminv(t[i], t[i + 1], t[i - 1], x[i], x[i + 1], x[i - 1])

    t_crit.append(_t)
    x_crit.append(_x)

t_crit = np.array(t_crit)
x_crit = np.array(x_crit)

# Alínea B
print("\n-> Alínea B")
amplitude = (F0 / MASSA) / (np.sqrt((OMEGA_F**2 - OMEGA_0**2)**2 + (B*OMEGA_F/MASSA)**2))
print(f"Amplitude: {amplitude} metros, quando wf = {OMEGA_F} rad/s")

# Alínea C
print("\n-> Alínea C")
OMEGA_F = 0.5
amplitude = (F0 / MASSA) / (np.sqrt((OMEGA_F**2 - OMEGA_0**2)**2 + (B*OMEGA_F/MASSA)**2))
print(f"Amplitude: {amplitude} metros, quando wf = {OMEGA_F} rad/s")
print(f"Neste caso temos que w0 = {OMEGA_0} rad/s é igual a wf = {OMEGA_F} rad/s, ou seja, na realidade estamos a calcular a amplitude quando tanto o w0 como o wf têm o mesmo valor. Quando isso acontece é porque estamos perante a ressonância, que é o valor máximo que a amplitude pode ter e acontece exatamente quando wf = w0")
print(f"Assim, {amplitude} metros não é só a amplitude quando wf = 0.5 rad/s mas também é o valor da ressonância, por isso é que a diferença entre os valores é tão alta.")
print("Este facto pode também ser explicado pela fórmula da amplitude \n\n(F0 / MASSA) / (np.sqrt((OMEGA_F**2 - OMEGA_0**2)**2 + (B*OMEGA_F/MASSA)**2))\n\nque será máxima quando o denominador for mínimo, ou seja quando \n\n(OMEGA_F**2 - OMEGA_0**2)**2\n\nfor mínimo, que só acontece quando wf = w0, obtendo-se a ressonância.")

plt.plot(t, x, "r")
plt.plot(t_crit, x_crit, "ob")

plt.xlabel("t (s)")
plt.ylabel("x (m)")
plt.title("Posição")
plt.show()

plt.plot(t, v, "b")
plt.xlabel("t (s)")
plt.ylabel("v (m/s)")
plt.title("Velocidade")
plt.show()