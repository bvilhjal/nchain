# Figure
import pylab as pl
import numpy as np
from matplotlib import pyplot as plt
from pylab import  get_cmap

Matrix_counts = np.genfromtxt ('test_presence_absence_matrix.csv', delimiter=",")
Matrix_counts = Matrix_counts[~(Matrix_counts==0).all(1)]

print 'The shape of the matrix of counts is:', Matrix_counts.shape

# Saving the matrix
#np.savetxt("foo.csv", Matrix_counts, delimiter=",")

#separate core and accessory genes:
plt.pcolor(Matrix_counts, cmap=get_cmap("RdBu"))
plt.ylim([0,199])
plt.title('Core and accessory genes of 200 strains')
plt.xlim([0,18760])
plt.xlabel('Genes')
plt.ylabel('Strains')
#plt.show()
plt.savefig('presence_absence_genospecies')
#plt.show()
