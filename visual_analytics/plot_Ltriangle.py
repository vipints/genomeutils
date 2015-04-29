
import matplotlib as mp
import pylab as pl
import pandas as pd 
import numpy as np

dis = pd.read_csv("test.csv", header=0)
C = np.tril(dis)
sim = 1-dis
C = np.tril(sim)
N = sim.shape[1]
C = np.ma.masked_array(C, C == 0)

A = np.array([(y, x) for x in range(N, -1, -1) for y in range(N + 1)])
t = np.array([[0.5, 1], [0.5, -1]])
A = np.dot(A, t)
X = A[:, 1].reshape(N + 1, N + 1)
Y = A[:, 0].reshape(N + 1, N + 1)
fig = pl.figure(figsize=(20,20))
ax = fig.add_subplot(121, frame_on=False, aspect=2.0)
ax.set_xticks([])
ax.set_yticks([])
caxes = pl.pcolormesh(X, Y, np.flipud(C), axes=ax)
ax.set_xlim(right=0)
outfile="similarity.pdf"
fig.savefig(outfile, bbox_inches='tight')
