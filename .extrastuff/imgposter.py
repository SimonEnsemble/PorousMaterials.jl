from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.gca(projection='3d')
data = np.load("data/db/SBMOF-1_coryimgposter.npy")
print(data.shape)
xl, yl, zl = data.shape

for i in range(5):
    for j in range(5):
        for k in range(5):
            if data[i,j,k] == 0:
                ax.text(i/5,j/5,k/5, str(data[i,j,k]), color = 'red')
            else:
                ax.text(i/5,j/5,k/5, str(data[i,j,k]), color = 'green')
#plt.axis('off')
ax.set_axis_off()
plt.savefig("imgposterasdf.png")
plt.show()
