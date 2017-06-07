# Figure 3A 
import pylab as pl
import numpy as np
import scipy
import matplotlib
from matplotlib import pyplot as plt
from pylab import  get_cmap
import seaborn as sns
import pandas as pd

#Matrix_counts = np.genfromtxt("presence_absence_matrix_zeros_ones.csv", delimiter=",")

import sys
sys.setrecursionlimit(100000)

def black_white_heatmap(Matrix_counts):
	
	# Trying clustermap 
	gsA = ['#386CB0']
	gsB = ['#FB8072']
	gsC = ['#1B9E77']
	gsD = ['#984EA3']
	gsE = ['#F0027F']

	colors = 33*gsA + 34*gsB + 116*gsC + 5*gsD + 11*gsE

	#col_linkage = scipy.cluster.hierarchy.linkage(Matrix_counts[0:199,0:14000], method = 'average', )
	
	pandas_data_frame = pd.DataFrame(Matrix_counts[0:199,0:14000])
	
	fig = plt.figure()
	g = sns.clustermap(pandas_data_frame, row_cluster=True, col_linkage = col_linkage)
	plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation= 0)
	fig.savefig('presence_absence_genospecies_zeros_ones_test.png')

#black_white_heatmap(Matrix_counts)

def len_forced(obj):
	ob_len = len(obj)
	if ob_len > 5:
		ob_len = 5
	return(ob_len)

def colorful_heatmap(Matrix_counts):
	
	# Removing one empty row
	Matrix_counts = Matrix_counts[~(Matrix_counts==0).all(1)]
	b = np.sum(Matrix_counts, axis = 0)
	idx = b.argsort()

	idx = sorted(range(0,18760), key=lambda x: list(np.unique(Matrix_counts[:,x])))
	idx = sorted(idx, key=lambda x: len_forced(list(np.unique(Matrix_counts[:,x]))))
	Matrix_counts = Matrix_counts[:,idx]

	print 'The shape of the matrix of counts is:', Matrix_counts.shape

	#countries = {'gsA': 1, 'gsB':2, 'gsC':3, 'gsD':4, 'gsE': 5}

	# Trying clustermap 
	gsA = ['#386CB0']
	gsB = ['#FB8072']
	gsC = ['#1B9E77']
	gsD = ['#984EA3']
	gsE = ['#F0027F']

	colors = 34*gsA + 34*gsB + 115*gsC + 5*gsD + 11*gsE

	#method='average', metric= 'hamming'
	pandas_data_frame = pd.DataFrame(Matrix_counts[0:199, 0:14000])
	g = sns.clustermap(pandas_data_frame, row_cluster=True, row_colors = colors)
	plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
	plt.savefig('presence_absence_genospecies_14000_rowcluster', format = 'pdf')
	
	#separate core and accessory genes:
	#plt.pcolor(Matrix_counts, cmap=matplotlib.colors.ListedColormap(colors))
	#plt.ylim([0,199])
	#plt.title('Core and accessory genes of 200 strains')
	#plt.xlim([0,18760])
	#plt.xlabel('Genes')
	#plt.ylabel('Strains')
	#plt.savefig('presence_absence_genospecies_white')
	#plt.show()
#Matrix_counts = np.genfromtxt("presence_absence_matrix.csv", delimiter=",")
#colorful_heatmap(Matrix_counts = Matrix_counts)