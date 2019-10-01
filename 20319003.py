import numpy as np
import matplotlib.pyplot as plt
#from matplotlib import gridspec

mu_a = -45.02
sigma_a = 6.01
mu_b = -31.79
sigma_b = 1.38

def distribution(x, mu, sigma):
    lebar = 1/np.power((2*np.pi*np.power((sigma),2)),0.5)
    tinggi = np.exp(-np.power((x-mu),2)/(2*np.power((sigma),2)))
    dist = lebar * tinggi
    return(dist)
    
xpoints_a = np.linspace(mu_a - 3*sigma_a, mu_a + 3*sigma_a, 100)
ypoints_a = distribution(xpoints_a, mu_a, sigma_a)
xpoints_b = np.linspace(mu_b - 3*sigma_b, mu_b + 3*sigma_a, 100)
ypoints_b = distribution(xpoints_b, mu_b, sigma_b)

#npoints = [100, 200, 300, 500, 700, 1000, 3000, 5000, 10000, 15000, 20000, 30000, 40000, 50000]
npoints = [1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000]

nbin = 30
b1 = np.linspace(-50, -40, nbin)
b2 = np.linspace(-38, -28, nbin)

biru1 = []
biru2 = []
for i in range (len(b1)):
    biru1.append(distribution(b1[i], mu_a, sigma_a))
    biru2.append(distribution(b2[i], mu_b, sigma_b))

u1 = [[] for i in range (len(npoints)-1)]
u2 = [[] for i in range (len(npoints)-1)]
z1 = [[] for i in range (len(npoints)-1)]
z2 = [[] for i in range (len(npoints)-1)]

n1 = [[] for i in range (len(npoints)-1)]
t1 = [[] for i in range (len(npoints)-1)]
luar1 = [[] for i in range (len(npoints)-1)]
jumlah1 = []
persen1 = []

n2 = [[] for i in range (len(npoints)-1)]
t2 = [[] for i in range (len(npoints)-1)]
luar2 = [[] for i in range (len(npoints)-1)]
jumlah2 = []
persen2 = []

for i in range(len(npoints)-1):
    u1[i] = np.random.random(npoints[i])
    u2[i] = np.random.random(npoints[i])

    z1[i] = 6 * np.sqrt(-2.0*np.log(u1[i])) * np.cos(2*np.pi*u2[i]) + mu_a
    z2[i] = np.sqrt(-2.0*np.log(u1[i])) * np.sin(2*np.pi*u2[i]) + mu_b
    
    a = np.histogram(z1[i], b1, density=True)
    n1[i] = a[0]
    hahaha = np.histogram(z2[i], b2, density=True)
    n2[i] = hahaha[0]
    
    for j in range(nbin-1):
        if n1[i][j] > biru1[j]:
            huft = n1[i][j] - biru1[j]
            luar1[i].append(huft)
        if n2[i][j] > biru2[j]:
            huftisasi = n2[i][j] - biru2[j]
            luar2[i].append(huftisasi)
    jumlah1.append(np.sum(luar1[i]))
    t1[i] = np.sum(n1[i])
    persen1.append(round(((t1[i] - jumlah1[i])/t1[i]), 3))
    jumlah2.append(np.sum(luar2[i]))
    t2[i] = np.sum(n2[i])
    persen2.append(round(((t2[i] - jumlah2[i])/t2[i]), 3))

for i in range(len(npoints)-1):
    fig = plt.figure(figsize=(6, 6))
    grid = plt.GridSpec(4, 4, hspace=0.1, wspace=0.1)
    main_ax = fig.add_subplot(grid[1:, :-1])
    y_hist = fig.add_subplot(grid[1:, -1], sharey=main_ax)
    x_hist = fig.add_subplot(grid[0, :-1], sharex=main_ax)
    
    # scatter points on the main axes
    main_ax.plot(z1[i], z2[i], 'ok', markersize=3, alpha=0.5)
    main_ax.set_ylim(-40, -25)
    main_ax.text(-60, -39, 'N = ' + str(npoints[i]), fontsize = 14)
    main_ax.set_xlabel(r'$z_1$, RV (LAMOST) in $kms^{-1}$')
    main_ax.set_ylabel(r'$z_2$, RV (Gaia) in $kms^{-1}$')
    
    # histogram on the attached axes
    x_hist.hist(z1[i], b1, histtype='stepfilled', orientation='vertical', color='gray', density=True)
    x_hist.plot(xpoints_a, ypoints_a)
    x_hist.set_ylim(0, 0.2)
    x_hist.text(-60, 0.1, 'ratio = ' + str(persen1[i]), fontsize = 13)
    x_hist.set_ylabel(r'$N_{normalized}$')
    plt.setp(x_hist.get_xticklabels(), visible=False)
    #x_hist.invert_yaxis()
    
    y_hist.hist(z2[i], b2, histtype='stepfilled', orientation='horizontal', color='gray', density=True)
    y_hist.plot(ypoints_b, xpoints_b)
    y_hist.set_xlim(0, 0.5)
    y_hist.text(0.3, -26, 'ratio = ' + str(persen2[i]), fontsize = 13, rotation=270)
    y_hist.set_xlabel(r'$N_{normalized}$')
    plt.setp(y_hist.get_yticklabels(), visible=False)
    #y_hist.invert_xaxis()
    #plt.show()
    filename = 'plots1000_' + str(i) +'.png'
    plt.savefig(filename)
    plt.close()
'''

fig = plt.figure(figsize=(6, 6))
grid = plt.GridSpec(4, 4, hspace=0.1, wspace=0.1)
main_ax = fig.add_subplot(grid[1:, :-1])
y_hist = fig.add_subplot(grid[1:, -1], sharey=main_ax)
x_hist = fig.add_subplot(grid[0, :-1], sharex=main_ax)
    
# scatter points on the main axes
main_ax.plot(z2[0], z1[0], 'ok', markersize=3, alpha=0.5)
main_ax.set_ylim(-40, -25)
#main_ax.text(-60, -39, 'N = ' + str(npoints[0]), fontsize = 14)
main_ax.set_xlabel(r'$z_2$, RV (LAMOST) in $kms^{-1}$')
main_ax.set_ylabel(r'$z_1$, RV (Gaia) in $kms^{-1}$')
    
# histogram on the attached axes
x_hist.hist(z1[0], b1, histtype='stepfilled', orientation='vertical', color='gray', density=True)
x_hist.plot(xpoints_b, ypoints_b)
x_hist.set_ylim(0, 0.5)
#x_hist.text(-60, 0.3, 'ratio = ' + str(persen2[0]), fontsize = 13)
x_hist.set_ylabel(r'$N_{normalized}$')
plt.setp(x_hist.get_xticklabels(), visible=False)
#x_hist.invert_yaxis()
    
y_hist.hist(z2[0], b2, histtype='stepfilled', orientation='horizontal', color='gray', density=True)
y_hist.plot(ypoints_a, xpoints_a)
y_hist.set_xlim(0, 0.5)
#y_hist.text(0.3, -26, 'ratio = ' + str(persen2[i]), fontsize = 13, rotation=270)
y_hist.set_xlabel(r'$N_{normalized}$')
plt.setp(y_hist.get_yticklabels(), visible=False)
#y_hist.invert_xaxis()
plt.show()
#filename = 'plots1000_' + str(i) +'.png'
#plt.savefig(filename)
#plt.close()'''