import numpy as np
import matplotlib.pyplot as plt

# Constantes

p = 350 # watts
x0 = 0
v0 = 1.0 # m/s
t0 = 0
m = 70 # kg
c_res = 0.9
u = 0.004
area = 0.35 # m^2
ro = 1.225
g = 9.81
tf = 100 # INSTANTE FINAL
dt = 0.001 # DELTA T

# Resolução

def EulerMethod (g, x0, v0, t0, tf, dt, theta):
    #variáveis
    n = int(np.ceil(tf/dt)) #número de passos (tempo total a dividir pelo tempo de cada passo)
    x = np.zeros(n + 1)
    v = np.zeros(n + 1)
    a = np.zeros(n + 1)
    t = np.zeros(n + 1)
    #inicialização
    x[0] = x0
    v[0] = v0
    t[0] = t0
    #loop de integração
    for i in range(n):
        a[i] = p/(m*v[i])-g*np.sin(theta)-(c_res*area*ro*v[i]*v[i]*np.cos(theta))/(2*m)-u*g*np.cos(theta)
        x[i+1] = x[i] + v[i] * dt
        v[i+1] = v[i] + a[i] * dt
        t[i+1] = t[i] + dt
    return n, x, v, a, t

n, x, v, a, t = EulerMethod (g, x0, v0, t0, tf, dt, np.deg2rad(0))

# Alínea A
print("\n-> Alínea A")
print(f"A velocidade terminal do ciclista é {v[-1]} m/s (ajustado para {t[-1]} segundos, podemos mudar o instante final na variável tf)")

# Alínea B
print("\n-> Alínea B")
for i in range(len(x)):
    if x[i] >= 500.0:
        print(f"O ciclista percorre {x[i]} metros ao final de {t[i]} segundos")
        print(f"Para além disso, tem uma velocidade de {v[i]} m/s neste instante")
        break

# Alínea C
print("\n-> Alínea C")
for i in range(len(x)):
    if x[i] >= 200.0:
        print(f"Quando atinge {x[i]} metros, o ciclista atinge uma velocidade de {v[i]} m/s ao final de {t[i]} segundos")
        print(f"A velocidade neste instante ({v[i]} m/s) representa {(v[i]/v[-1])*100}% da velocidade terminal ({v[-1]} m/s))")
        break

print(f"\nA aceleração em durante a velocidade terminal é igual a {a[-1]}, como seria de esperar")

# Plot
plt.title("Evolução temporal da posição do ciclista")
plt.plot(t,x, label="Posição")
plt.xlabel("Tempo (s)")
plt.ylabel("Posição (m)")
plt.legend()
plt.show()

plt.title("Evolução temporal da velocidade do ciclista")
plt.plot(t,v, label="Velocidade")
plt.xlabel("Tempo (s)")
plt.ylabel("Velocidade (m/s)")
plt.legend()
plt.show()

plt.title("Evolução temporal da aceleração do ciclista")
plt.plot(t,a, label="Aceleração")
plt.xlabel("Tempo (s)")
plt.ylabel("Aceleração (m/s^2)")
plt.legend()
plt.show()