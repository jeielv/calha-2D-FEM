import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.constants import epsilon_0

def auto_fmt(x, casas=1, lim=1e-3):
	if x == 0:
		return "0.0"
	if abs(x) >= lim and abs(x) < 1/lim:
		return f"{x:.{casas}f}"
	else:
		return f"{x:.{casas}e}"

def q(i, elem):
	if i == 0:
		return (elem[1].y - elem[2].y)
	elif i == 1:
		return (elem[2].y - elem[0].y)
	else: return (elem[0].y - elem[1].y)

def r(i, elem):
	if i == 0:
		return (elem[2].x - elem[1].x)
	elif i == 1:
		return (elem[0].x - elem[2].x)
	else: return (elem[1].x - elem[0].x)

def K(i, j, elem, D, e):
	#calculo K da forma convencional, usando os q e r, tambem calculados convencionalmente nas funções definidas acima
	qi = q(i, elem)
	ri = r(i, elem)
	qj = q(j, elem)
	rj = r(j, elem)

	K = (e / (2 * D)) * (qi * qj + ri * rj)

	return K

def calc_sist_global(elementos, Nx, Ny):

	#o valor D (o dobro da area do elemento triangular) eh o mesmo para todos os elementos
	#portanto, pega-se somente o primeiro elemento msm:
	elemento = elementos[0]
	D = np.linalg.det([[1, elemento[0].x, elemento[0].y], [1, elemento[1].x, elemento[1].y], [1, elemento[2].x, elemento[2].y]])

	#abaixo, utiliza-se a biblioteca scipy para pegar a permissividade do ar, que sera utilizada para calcular K
	#de todos os elementos (todos sao preenchidos por ar)
	e = epsilon_0

	sist_glob = []
	#inicializo a matriz pondo tudo 0
	for i in range(Nx * Ny):
		linha = []
		for j in range (Nx * Ny):
			linha.append(0.0)
		sist_glob.append(linha)

	#no codigo abaixo, calculo cada K da matriz elementar de cada elemento 
	#(identificado pelos seus nós) e ja incluo esses valores K no sistema global

	#a inclusao dos valores no sistema global foi feita com o seguinte pensamento: a matriz global
	#comeca toda zerada e vai se somando cumulativamente às posições os valores vindos das matrizes
	#elementares. dado um valor Kij em uma matriz elementar associada ao elemento k, esse valor sera
	#incluido no sistema global como uma soma à posição cuja linha tem indice igual ao indice do nó
	#global correspondente ao i-esimo nó local do elemento k e cuja coluna tem indice igual ao indice do nó
	#global correspondente ao j-esimo nó local do elemento k
	for elem in elementos:
		for i in range(3):
			for j in range(3):
				Kij = K(i, j, elem, D, e)
				sist_glob[elem[i].idx][elem[j].idx] += Kij

	#aqui, altero as linhas do sistema global correspondente aos nós das paredes da calha
	#de forma a respeitar as condicoes de contorno
	# i == cnt * Nx: nós na parede esquerda da calha
	# i < Nx: nós na parede debaixo da calha
	# i == cnt * Nx - 1: nós na parede direita da calha
	# cnt == Ny: nós na parede de cima da calha
	# todos esses nós vao ter suas linhas do sistema global zeradas, exceto qnd i == j, que valera 1.0
	cnt = 0
	for i in range(Nx * Ny):
		linha = []
		if i == cnt * Nx or i < Nx or i == cnt * Nx - 1 or cnt == Ny:
			for j in range (Nx * Ny):
				sist_glob[i][j] = (0.0 if i != j else 1.0)
			if i == cnt * Nx: cnt += 1
	
	return sist_glob

def melhor_nx_ny(N):
	# razao do dominio (L / (L/2) = 2)
	r = 2.0

	# M deve ser N / 2, que eh a qtd de triangulos retangulos; se N impar, ent arredondo
	if N % 2 != 0:
		M = int(N // 2)
	else:
		M = N // 2

	melhor_erro = 999999
	melhor_s = 0
	melhor_t = 0

	#procuro todos os pares s, t com s * t = M (s = Nx - 1 e t = Ny - 1) e escolho o par com menor erro
	#usei uma forma de evitar que Nx e Ny saiam muito desbalanceados quando N / 2 eh primo
	for s in range(1, int(math.sqrt(M))*4 + 5):
			for t in range(1, int(math.sqrt(M))*4 + 5):

				erro_geo = abs(s/t - r)
				erro_area = abs(s * t - M)

				erro_total = 2 * erro_geo + erro_area

				if erro_total < melhor_erro:
					melhor_erro = erro_total
					melhor_s = s
					melhor_t = t

	Nx = melhor_s + 1
	Ny = melhor_t + 1

	#com esse algoritmo, pode haver casos em que serão dados valores bem distribuídos de Nx e Ny, mas com um valor total 
	#de elementos nao exatamente igual a N, mas próximo disso
	return Nx, Ny

class Nos:
	def __init__(self, x, y, v, idx):
		self.x = x
		self.y = y
		self.v = v
		self.idx = idx

N = int(input("quantos elementos? "))
V0 = float(input("de um valor para V0 (em V): "))
L = float(input("de um valor para L (em metros): "))

#Nx = qtd de nós no eixo X; Ny = qtd de nós no eixo Y
#usei um algoritmo para encontrar valores de Nx e Ny que nn deixem mt espaco sem nó nos eixos

Nx, Ny = melhor_nx_ny(N)

#valores de variacao de x e de y para cada nó
dx = L / (Nx - 1)
dy = L / (2 * (Ny - 1))

#aqui, faco a matriz que vai conter todos os nós e suas informacoes, simulando o plano cartesiano
matriz_nos = []
for i in range(Ny):
	aux = []
	for j in range(Nx):
		#perceba que preencho inicialmente quase todos os nós com V = 0, exceto os da parede esquerda da calha
		aux.append(Nos(j * dx, i * dy, 0 if j != 0 else V0, Nx * i + j))
	matriz_nos.append(aux)

elementos = []

#aq, tenho a qtd real de elementos (pois dependendo do N dado, foi aproximado)
qtd_el = 2 * (Nx - 1) * (Ny - 1)

#abaixo, crio a matriz onde cada linha representa um elemento
for i in range(Ny - 1):
	for j in range(Nx - 1):
		aux = [matriz_nos[i][j], matriz_nos[i][j + 1], matriz_nos[i + 1][j]]
		elementos.append(aux)
		aux = [matriz_nos[i][j + 1], matriz_nos[i + 1][j + 1], matriz_nos[i + 1][j]]
		elementos.append(aux)

#aqui, calculo o sistema global K tal que K * V = F (V eh a matriz com as variaves q representam as tensoes nos nós)
sist_global = calc_sist_global(elementos, Nx, Ny)

#como nao ha densidade volumetrica no dominio, F tem quase todos os valores 0, exceto na linha referente aos valores de contorno
F = []
cnt = 0
for i in range(Nx * Ny):
	aux = []
	if i == cnt * Nx:
		aux = [V0]
		cnt += 1
	else:
		aux = [0.0]
	F.append(aux)

#na linha abaixo, resolvo o sistema para encontrar a matriz coluna contendo os valores de V para cada nó
Vc = np.linalg.solve(sist_global, F)

#aqui, so por questao de organizacao, preencho os valores certos dos potenciais na matriz de nós
for i in range(Ny):
	for j in range(Nx):
		matriz_nos[i][j].v = Vc[i * Nx + j][0]

#nas linhas abaixo, crio as matrizes X, Y e V que serão utilizadas para plotar os gráficos
X = []
Y = []
V = []

for i in range(Ny):
	auxX = []
	auxY = []
	auxV = []
	for j in range(Nx):
		auxX.append(matriz_nos[i][j].x)
		auxY.append(matriz_nos[i][j].y)
		auxV.append(matriz_nos[i][j].v)
	X.append(auxX)
	Y.append(auxY)
	V.append(auxV)

X = np.array(X)
Y = np.array(Y)
V = np.array(V)

########## DISTRIBUICAO DO POTENCIAL EM GRÁFICO 3D ##########

fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')

ax.plot_surface(X, Y, V)

plt.xlabel("x (m)")
plt.ylabel("y (m)")
ax.set_zlabel("Potencial V")

plt.show()

#############################################################

######## DISTRIBUICAO DO POTENCIAL EM GRAFICO DE COR ########

plt.figure()
plt.pcolormesh(X, Y, V, shading='auto')
plt.colorbar(label="Potencial V")

plt.xlabel("x (m)")
plt.ylabel("y (m)")

plt.title("Distribuição do Potencial na Calha")
plt.show()

##############################################################

########## DISTRIBUICAO DO POTENCIAL EM GRÁFICO MANTENDO X FIXO ###########

j = Nx // 2

plt.plot(Y[:, j], V[:, j])
plt.xlabel("y (m)")
plt.ylabel("Potencial V (V)")
plt.title(f"Variação do potencial ao longo de y para x = {X[0][j]:.3f} m")
plt.grid()
plt.show()

###########################################################################

########## DISTRIBUICAO DO POTENCIAL EM GRÁFICO MANTENDO Y FIXO ##########

i = Ny // 2

plt.plot(X[i], V[i])
plt.xlabel("x (m)")
plt.ylabel("Potencial V (V)")
plt.title(f"Variação do potencial ao longo de x para y = {Y[i][0]:.3f} m")
plt.grid()
plt.show()

###########################################################################