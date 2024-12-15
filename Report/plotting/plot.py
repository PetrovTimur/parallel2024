import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
from matplotlib.ticker import ScalarFormatter, MaxNLocator

ops = ['Dot', 'AXpY', 'SpMV', 'CGSolve']

def plotSingle(x, y, xlabel=None, ylabel=None, title=''):
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(8, 8), sharey=True)
    colors='cmyk'
    for i in range(4):
        axs.plot(x, y[i], marker='o', markeredgecolor=colors[i], c=colors[i])
    axs.set_xscale('log')
    axs.set_xlabel(xlabel)
    axs.set_ylabel(ylabel)
    axs.set_title(title)
    axs.legend(ops)
    # axs.xaxis.set_major_formatter(ScalarFormatter(useOffset=True))
    axs.set_xticks(x)

    plt.show()


def plotMulti(x, y, xlabel=None, ylabel=None, title=''):
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(8, 8), sharey=True)
    colors = 'cmyk'
    for i in range(3):
        axs.plot(x, y[i], marker='.', markeredgecolor=colors[i], c=colors[i])
        # axs[i].set_xscale('log')
    # axs.set_xscale('log', base=2)
    axs.set_xlabel(xlabel)
    axs.set_ylabel(ylabel)
    # axs.set_title(title)
    axs.legend([f'1e{i}' for i in range(5, 9)])
    axs.xaxis.set_major_formatter(ScalarFormatter())
    axs.set_xticks(x)

    # axs.plot([2, 32], [y[-1][0], 16 * y[-1][0]])
    plt.show()


def plotSingle2(x, y, xlabel=None, ylabel=None, title=''):
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(8, 8), sharey=True)
    colors='cmyk'
    for i in range(4):
        axs.plot(x, y[i], marker='o', markeredgecolor=colors[i], c=colors[i])
    axs.set_xscale('log')
    axs.set_yscale('log')
    axs.set_xlabel(xlabel)
    axs.set_ylabel(ylabel)
    axs.set_title(title)
    axs.legend(ops)
    # axs.xaxis.set_major_formatter(ScalarFormatter(useOffset=True))
    axs.set_xticks(x)

    plt.show()


# data_init = np.array([[0.0593148, 0.0322724, 0.0157722, 0.0112697, 0.0129851, 0.00893087]])
# plotMulti([1, 2, 4, 6, 8, 12], data_init[0][0] / data_init, xlabel='Threads', ylabel='Ускорение')
#
# data_init = np.array([[1.09502, 0.566558, 0.29257, 0.194787, 0.190751, 0.131433]])
# plotMulti([1, 2, 4, 6, 8, 12], data_init[0][0] / data_init, xlabel='Threads', ylabel='Ускорение')

# single_thread_memory = np.array([[1.8, 16.3, 160.2, 1568.8], [2.6, 24.2, 240.3, 2350.1], [2.2, 19.2, 197.3, 1929.8], [3.1, 28.1, 277.3, 2711.0]])
# plotSingle2([1e5, 1e6, 1e7, 1e8], single_thread_memory, xlabel='Размер данных', ylabel='Память, MB')

# single_thread_time = np.array([[0.00198042, 0.022299, 0.156074, 1.86033], [0.00181884, 0.0291459, 0.359465, 3.66604], [0.0051881, 0.0505653, 0.506143, 5.06442], [0.117945, 1.16356, 10.3548, 109.142]])
# plotSingle2([1e5, 1e6, 1e7, 1e8], single_thread_time, xlabel='Размер данных', ylabel='Время, s')

# single_data_gflops = np.array([[0.971881, 0.971756, 0.964252, 0.961759], [0.63664, 0.614648, 0.630975, 0.618567], [0.253414, 0.258222, 0.262578, 0.262851], [0.190021, 0.197509, 0.204131, 0.193596]])
# plotSingle([1e5, 1e6, 1e7, 1e8], single_data_gflops, xlabel='Размер данных', ylabel='GFLOPS')

# multithread_dot_data_gflops = np.array([[1.66884, 1.72424, 1.70368, 1.74324], [3.01223, 3.33049, 3.35082, 3.31488], [4.1214, 4.9362, 5.06217, 4.94426],
#                                         [3.59558, 3.78659, 3.8511, 3.88338], [3.98398, 4.16773, 4.24283, 4.25796]]).T
# plotMulti([2, 4, 6, 8, 12], multithread_dot_data_gflops, xlabel='Threads', ylabel='GFLOPS', title='Dot op')
#
# multithread_axpy_data_gflops = np.array([[1.01575, 0.981598, 1.01331, 0.973743], [1.85591, 1.80562, 1.64489, 2.04797], [2.92658, 3.25297, 2.96159, 3.45678],
#                                          [2.63385, 2.66633, 2.74976, 3.0151], [2.80172, 3.05359, 3.20047, 3.24413]]).T
# plotMulti([2, 4, 6, 8, 12], multithread_axpy_data_gflops, xlabel='Threads', ylabel='GFLOPS', title='AXpY op')
#
# multithread_spmv_data_gflops = np.array([[0.44424, 0.505814, 0.508037, 0.510398], [0.896403, 0.985598, 0.974944, 1.0025,], [1.39135, 1.48301, 1.50278, 1.50783],
#                                          [1.48747, 1.49436, 1.49074, 1.53363], [1.86673, 2.09183, 2.20504, 2.25376,]]).T
# plotMulti([2, 4, 6, 8, 12], multithread_spmv_data_gflops, xlabel='Threads', ylabel='GFLOPS', title='SpMV op')
#
# multithread_cgsolve_data_gflops = np.array([[0.353369, 0.361824, 0.374478, 0.364689], [0.642523, 0.697794, 0.713056, 0.722238], [0.926049, 1.01528, 0.999607, 1.04626],
#                                             [0.784544, 0.839341, 0.876031, 0.878439], [1.07405, 1.13838, 1.24984, 1.27949]]).T
# plotMulti([2, 4, 6, 8, 12], multithread_cgsolve_data_gflops, xlabel='Threads', ylabel='GFLOPS', title='CGSolve op')


# np.save('desktop', np.array([single_data_gflops, multithread_spmv_data_gflops.T, multithread_axpy_data_gflops.T, multithread_spmv_data_gflops.T]))

# ---------------------------------------------------------------------------POLUS---------------------------------------------------------------------

# single_data_gflops = np.array([[0.151781, 0.151825, 0.151764, 0.149401], [0.111303, 0.1113, 0.111051, 0.105555], [0.0639949, 0.0639803, 0.0638165, 0.063684], [0.0418468, 0.041932, 0.0419215, 0.0419531]])
# # plotSingle([1e5, 1e6, 1e7, 1e8], single_data_gflops, xlabel='Размер данных', ylabel='GFLOPS')
# #
# multithread_dot_data_gflops = np.array([[0.302481, 0.301435, 0.301008, 0.297214], [0.591624, 0.596378, 0.596087, 0.586329], [1.06538, 1.10087, 1.13066, 1.11105],
#                                         [1.50253, 1.65331, 1.67899, 1.55883], [1.76894, 2.18126, 2.21464, 2.14167], [2.09552, 2.64876, 2.71706, 2.60546],
#                                         [1.86439, 2.7102, 2.85122, 2.83682], [1.60101, 2.74435, 3.00793, 2.98526]]).T
# plotMulti([2, 4, 8, 16, 32, 40, 60, 80], multithread_dot_data_gflops, xlabel='Threads', ylabel='GFLOPS', title='Dot op')
#
# multithread_axpy_data_gflops = np.array([[0.226541, 0.225149, 0.225855, 0.214907], [0.435156, 0.437415, 0.436452, 0.410406], [0.820984, 0.829154, 0.829098, 0.780887],
#                                          [1.18123, 1.24651, 1.25766, 1.13494], [1.4761, 1.64346, 1.65129, 1.59584], [1.78809, 2.00608, 2.01414, 1.92499],
#                                          [1.79211, 2.05768, 2.10354, 2.08548], [1.49612, 2.0882, 2.17034, 2.16871,]]).T
# plotMulti([2, 4, 8, 16, 32, 40, 60, 80], multithread_axpy_data_gflops, xlabel='Threads', ylabel='GFLOPS', title='AXpY op')
# #
# multithread_spmv_data_gflops = np.array([[0.127937, 0.128035, 0.128091, 0.127949], [0.250717, 0.251794, 0.251906, 0.251918], [0.479233, 0.478283, 0.48203, 0.48181],
#                                          [0.719244, 0.728821, 0.728144, 0.727078], [0.973038, 0.982708, 0.983604, 0.983612], [1.18239, 1.19956, 1.19805, 1.19642],
#                                          [1.23793, 1.26545, 1.28508, 1.28806,], [1.20777, 1.28013, 1.29833, 1.30325]]).T
# plotMulti([2, 4, 8, 16, 32, 40, 60, 80], multithread_spmv_data_gflops, xlabel='Threads', ylabel='GFLOPS', title='SpMV op')
#
# multithread_cgsolve_data_gflops = np.array([[0.0834881, 0.0834048, 0.0833627, 0.0833596 ], [0.158126, 0.163314, 0.163861, 0.16394], [0.292248, 0.306937, 0.308487, 0.308656],
#                                          [0.397054, 0.453136, 0.459019, 0.460665], [0.454436, 0.590692, 0.60869, 0.61403], [0.535829, 0.70771, 0.73203, 0.736887],
#                                          [0.686739, 0.732487, 0.772683, 0.795443], [0.359764, 0.702728, 0.787378, 0.804724]]).T
# plotMulti([2, 4, 8, 16, 32, 40, 60, 80], multithread_cgsolve_data_gflops, xlabel='Threads', ylabel='GFLOPS', title='CGSolve op')
#

# multithread_dot_data_gflops = np.array([[0.300312, 0.300922, 0.300709, 0.296554,], [0.598482, 0.601168, 0.602436, 0.591719], [1.07451, 1.12729, 1.12585, 1.1079,],
#                                         [1.81628, 2.16618, 2.17559, 2.1078,], [2.03445, 3.10261, 3.2463, 3.00299,], [2.41637, 3.81488, 4.08225, 3.80341,],
#                                         [1.95825, 3.82272, 4.42505, 4.34894, ], [1.76719, 4.52494, 5.28819, 5.2221,], [1.35162, 3.71911, 4.99139, 5.60368,], [1.10429, 4.03472, 5.67538, 5.8043,]]).T
# plotMulti([2, 4, 8, 16, 32, 40, 60, 80, 120, 160], multithread_dot_data_gflops, xlabel='Threads', ylabel='GFLOPS', title='Dot op')
#
# multithread_axpy_data_gflops = np.array([[0.221842, 0.220587, 0.220487, 0.209703,], [0.435063, 0.436763, 0.436876, 0.411089,], [0.811505, 0.828372, 0.828339, 0.778485,],
#                                          [1.45428, 1.60079, 1.60646, 1.52092,], [1.85669, 2.39876, 2.42714, 2.19839, ], [2.19982, 2.96329, 3.05072, 2.76473,],
#                                          [2.16292, 3.05818, 3.30163, 3.21257, ], [1.82774, 3.6345, 3.94273, 3.79492,], [1.68805, 3.71407, 4.10764, 4.11398,], [1.22549, 3.29921, 4.13518, 4.24579,]]).T
# plotMulti([2, 4, 8, 16, 32, 40, 60, 80, 120, 160], multithread_axpy_data_gflops, xlabel='Threads', ylabel='GFLOPS', title='AXpY op')
# #
# multithread_spmv_data_gflops = np.array([[0.126828, 0.127846, 0.127593, 0.127003,], [0.245407, 0.251102, 0.249896, 0.250845,], [0.465388, 0.48134, 0.482144, 0.478431,],
#                                          [0.859945, 0.934143, 0.938284, 0.844518,], [1.21996, 1.40484, 1.4312, 1.41582,], [1.40376, 1.75595, 1.78706, 1.75265,],
#                                          [1.43944, 1.85922, 2.01231, 2.00478,], [1.47216, 2.26789, 2.26497, 2.34695,], [1.49166, 2.36065, 2.5283, 2.52529,], [1.11245, 2.25598, 2.54125, 2.58539,]]).T
# plotMulti([2, 4, 8, 16, 32, 40, 60, 80, 120, 160], multithread_spmv_data_gflops, xlabel='Threads', ylabel='GFLOPS', title='SpMV op')
#
# multithread_cgsolve_data_gflops = np.array([[0.0826346, 0.0829105, 0.0830233, 0.0829571,], [0.155955, 0.160717, 0.161569, 0.16192,], [0.28055, 0.297274, 0.300877, 0.301486,],
#                                          [0.463826, 0.535491, 0.560002, 0.563258,], [0.599644, 0.76496, 0.773002, 0.823514,], [0.57977, 0.890344, 0.983785, 0.895145,],
#                                          [0.505283, 0.898503, 1.05905, 1.09827,], [0.550955, 1.05815, 1.23718, 1.25988,], [0.347093, 0.976246, 1.27145, 1.33129,], [0.313439, 0.955972, 1.28518, 1.354,]]).T
# plotMulti([2, 4, 8, 16, 32, 40, 60, 80, 120, 160], multithread_cgsolve_data_gflops, xlabel='Threads', ylabel='GFLOPS', title='CGSolve op')


# dot_speedup = multithread_dot_data_gflops.T / single_data_gflops[0]
# data = dot_speedup.T[-1]
# print(*[f"{spdup:.2f}" for spdup in data], sep=' & ')
#
# axpy_speedup = multithread_axpy_data_gflops.T / single_data_gflops[1]
# data = axpy_speedup.T[-1]
# print(*[f"{spdup:.2f}" for spdup in data], sep=' & ')
#
# spmv_speedup = multithread_spmv_data_gflops.T / single_data_gflops[2]
# data = spmv_speedup.T[-1]
# print(*[f"{spdup:.2f}" for spdup in data], sep=' & ')
#
# cgsolve_speedup = multithread_cgsolve_data_gflops.T / single_data_gflops[3]
# data = cgsolve_speedup.T[-1]
# print(*[f"{spdup:.2f}" for spdup in data], sep=' & ')


import pandas as pd

df = pd.read_csv('mpi_dot.csv', index_col=False, sep=',')
# print(df[df['Op'] == 'Dot'][['1e5', '1e6']].values)
plotMulti(df[df['Op'] == 'Dot']['NumProc'], np.array(df[df['Op'] == 'Dot'][['1e6', '1e7', '1e8']].values).T, xlabel='numProc', ylabel='GFLOPS', title='Dot')
plotMulti(df[df['Op'] == 'AXpY']['NumProc'], np.array(df[df['Op'] == 'AXpY'][['1e6', '1e7', '1e8']].values).T, xlabel='numProc', ylabel='GFLOPS', title='AXpY')
plotMulti(df[df['Op'] == 'SpMV']['NumProc'], np.array(df[df['Op'] == 'SpMV'][['1e6', '1e7', '1e8']].values).T, xlabel='numProc', ylabel='GFLOPS', title='SpMV')
plotMulti(df[df['Op'] == 'CGSolver']['NumProc'], np.array(df[df['Op'] == 'CGSolver'][['1e6', '1e7', '1e8']].values).T, xlabel='numProc', ylabel='GFLOPS', title='CGSolver')
# plotMulti([2, 4, 8, 16, 32, 40, 60, 80], multithread_dot_data_gflops, xlabel='Threads', ylabel='GFLOPS', title='Dot op')

