import numpy as N


nk = 3
nw = 12
k =  complex(1.,0.)*N.random.random(nk)+complex(0.,1.)*N.random.random(nk)
w =  complex(1.,0.)*N.random.random(nw)+complex(0.,1.)*N.random.random(nw)

fac = complex(1.,0.)*N.random.random(nk)+complex(0.,1.)*N.random.random(nk)


K = k[:,N.newaxis]

z = w/(K**2-w**2)

result = fac[:,N.newaxis]*z


dumb_result = 1j*N.zeros([len(k),len(w)])


for i in N.arange(len(k)):
    for j in N.arange(len(w)):

        dumb_result[i,j] = w[j]/(k[i]**2-w[j]**2)*fac[i]

error = N.linalg.norm(dumb_result-result)
