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
        axs.plot(x[0]*x[1], y[i], marker='.', markeredgecolor=colors[i], c=colors[i])
        # axs[i].set_xscale('log')
    axs.set_xscale('log', base=2)
    axs.set_yscale('log', base=2)
    axs.set_xlabel(xlabel)
    plt.tick_params('x', labelsize=10)
    axs.set_ylabel(ylabel)
    axs.set_title(title)
    axs.legend([f'1e{i}' for i in range(6, 9)])
    axs.xaxis.set_major_formatter(ScalarFormatter())
    labels = [f'{proc}' for proc, thread in x.T]
    # print(labels)
    # axs.xaxis.set_major_formatter(FormatStrFormatter())
    axs.yaxis.set_major_formatter(ScalarFormatter())
    # axs.yaxis.set_major_locator(MaxNLocator(integer=True))
    axs.set_xticks(x[0]*x[1])
    axs.set_xticklabels(labels)

    # axs.axvline(x=40, ymin=0.02, ymax=0.98, linestyle='--', label='test', c='red')

    fig.tight_layout()
    plt.show()


def plotCompare(x, y, xlabel=None, ylabel=None, title='', labels=None):
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
    colors = 'rgb'
    for i in range(2):
        axs.plot(x[i], y[i], marker='.', markeredgecolor=colors[i], c=colors[i])

    # axs.set_xscale('log', base=2)
    # axs.set_yscale('log', base=2)

    axs.set_xlabel('Threads (OpenMP)')
    axs.set_ylabel(ylabel)
    axs.set_title(title)
    axs.legend(['OpenMP', 'MPI'])

    # ticks = set(x[0]) | set(x[1])
    # print(ticks, x[0], x[1])
    axs.set_xticks(x[0])
    axs.set_xticklabels(x[0])
    axs.xaxis.set_major_formatter(ScalarFormatter())
    # axs.yaxis.set_major_locator(MaxNLocator(integer=True))
    plt.tick_params('x', labelsize=6)

    secax = axs.secondary_xaxis('top')
    secax.set_xlabel('numProc x Threads (MPI)')
    secax.minorticks_off()
    secax.set_xlim(axs.get_xlim())
    secax.set_xticks(x[1])
    secax.set_xticklabels(x[1])
    secax.xaxis.set_major_formatter(ScalarFormatter())
    secax.tick_params('x', labelsize=6)

    labels = [f'{proc}x{thread}' for proc, thread in labels.T]
    secax.set_xticklabels(labels)

    labels1, labels2 = axs.xaxis.get_ticklabels(), secax.xaxis.get_ticklabels()
    for i in [0, 1, 2, 4]:
        labels1[i].set_visible(False)
        labels2[i].set_visible(False)

    fig.tight_layout()
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


# np.save('desktop', np.array([single_data_gflops, multithread_spmv_data_gflops.T, multithread_axpy_data_gflops.T, multithread_spmv_data_gflops.T]))

# ---------------------------------------------------------------------------POLUS---------------------------------------------------------------------

single_data_gflops = np.array([[0.151781, 0.151825, 0.151764, 0.149401], [0.111303, 0.1113, 0.111051, 0.105555], [0.0639949, 0.0639803, 0.0638165, 0.063684], [0.0418468, 0.041932, 0.0419215, 0.0419531]])


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

file_path = 'time2.csv'
df = pd.read_csv(file_path, index_col=False, sep=',')
df_mpi = df[df['Type'] == 'MPI']
df_omp = df[df['Type'] == 'OpenMP']
df_cuda = df[df['Type'] == 'CUDA']

# fig, axs = plt.subplots(4, 1, sharex=True, figsize=(8, 16))
# for i, op in enumerate(['Dot', 'AXpY', 'SpMV', 'CGSolver']):
    # plotMulti(np.stack((df_mpi[df_mpi['Op'] == op]['NumProc'].to_numpy(), df_mpi[df_mpi['Op'] == op]['Threads'].to_numpy())), np.array(df_mpi[df_mpi['Op'] == op][['1e6', '1e7', '1e8']].values).T, xlabel='numProc', ylabel='GFLOPS', title=op)
    # plotMulti(np.stack((df_omp[df_omp['Op'] == op]['NumProc'].to_numpy(), df_omp[df_omp['Op'] == op]['Threads'].to_numpy())), np.array(df_omp[df_omp['Op'] == op][['1e6', '1e7', '1e8']].values).T, xlabel='Threads', ylabel='GFLOPS', title=op)

# for i, op in enumerate(['Dot', 'AXpY', 'SpMV', 'CGSolver']):
#     # print(type((df_omp[df_omp['Op'] == op]['NumProc']*df_omp[df_omp['Op'] == op]['Threads']).values))
#     x = [(df_omp[df_omp['Op'] == op]['NumProc']*df_omp[df_omp['Op'] == op]['Threads']).values, (df_mpi[df_mpi['Op'] == op]['NumProc']*df_mpi[df_mpi['Op'] == op]['Threads']).values]
#     y = [np.array(df_omp[df_omp['Op'] == op][['1e8']].values).squeeze(1), np.array(df_mpi[df_mpi['Op'] == op][['1e8']].values).squeeze(1)]
#     # print(y)
#     # print(x[0].shape, y[0].shape)
#     labels = np.stack((df_mpi[df_mpi['Op'] == op]['NumProc'].to_numpy(), df_mpi[df_mpi['Op'] == op]['Threads'].to_numpy()))
#     plotCompare(x, y, xlabel='numProc', ylabel='GFLOPS', title=op, labels=labels)

# for i, op in enumerate(['Dot', 'AXpY', 'SpMV', 'CGSolver']):
    # dot_speedup = np.array(df_mpi[df_mpi['Op'] == op][['1e6', '1e7', '1e8']].values) / df_mpi[df_mpi['Op'] == op][['1e6', '1e7', '1e8']].values[0]
    # dot_speedup = np.array(df_omp[df_omp['Op'] == op][['1e6', '1e7', '1e8']].values) / df_omp[df_omp['Op'] == op][['1e6', '1e7', '1e8']].values[0]
    # data = dot_speedup.T[-1]
    # print(*[f"{spdup:.2f}" for spdup in data], sep=' & ')

def plot_gflops_vs_size(data, sizes, xlabel='Size', ylabel='GFLOPS', title=''):
    x = [int(float(s)) for s in sizes]
    plt.figure(figsize=(8, 8))
    plt.plot(x, data, marker='o')
    plt.xscale('log')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    # plt.grid(True)
    plt.ylim(0, 1.2 * max(data))
    plt.show()

# for i, op in enumerate(['Dot', 'AXpY', 'SpMV', 'CGSolver']):
#     plot_gflops_vs_size(df_cuda[df_cuda['Op'] == op][['1e6', '1e7', '1e8']].values.T, ['1e6', '1e7', '1e8'])

construct = np.array([2.01431, 1.69289, 1.347, 1.12794, 0.87568])
transpose = np.array([1.35265, 0.652406, 0.360206, 0.235345, 0.186395])
adjacency = np.array([31.4697, 16.836, 8.59597, 6.33709, 4.94728])
fill = np.array([12.1859, 6.33143, 3.27431, 2.13725, 1.40589])

stages = np.stack([construct, transpose, adjacency, fill])

# print(stages)

from matplotlib.ticker import ScalarFormatter
def plot_time_vs_threads(data, threads, xlabel='Threads', ylabel='Speedup', title=''):
    plt.figure(figsize=(8, 8))
    plt.plot(threads, data, marker='o')
    plt.xscale('log', base=2)
    plt.yscale('log')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    # plt.yticks(data, [f'{x:.2f}' for x in data])
    formatter = ScalarFormatter(useMathText=False, useOffset=False)
    formatter.set_scientific(False)
    plt.gca().yaxis.set_major_formatter(formatter)
    plt.gca().yaxis.set_minor_formatter(formatter)
    plt.xticks(threads, [str(int(t)) for t in threads])
    plt.show()

for i, op in enumerate(['EN', 'NE', 'EE', 'Fill']):
    plot_time_vs_threads(stages[i][0] / stages[i], [1, 2, 4, 6, 12])