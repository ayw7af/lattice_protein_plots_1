#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# =============================================================================
# Begins plot.py
#
# This module visualize the trajectory outputs using the matplotlib package.
# =============================================================================

def running_mean(x, N):
    import numpy as np
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)

def gyration_calculation(conf_list, protein_name, T):
    from src.gyration import gyration

    if protein_name == "mer48A":
        text = "48-mer-A"
        from proteins.mer48A import seq
    elif protein_name == "mer27":
        text = "27-mer"
        from proteins.mer27 import seq
    else:
        text = "chignolin"
        from proteins.chignolin import seq

    steps = len(conf_list)
    stats = int(steps/1000)
    temp_conf = conf_list[::stats]
    r_g = []
    for i in temp_conf:
        r_g.append(gyration(i, seq))

    import numpy as np
    import matplotlib.pyplot as plt
    fig9 = plt.figure()
    label9 = "Radius of gyration of " + text + " at $T$ = " + str(T) + " K"
    plt.scatter(np.arange(0.0, steps, stats), r_g,
                c='black', s=1, marker='.')
    plt.scatter(np.arange(50 * stats, steps, stats),
                running_mean(r_g, 51),
                c='blue', s=0.5, marker='.')
    plt.title(label9)
    plt.ylabel("$R_g$($\mathrm{\AA}$)")
    plt.xlabel("Monte Carlo time (steps)")
    #plt.xlim(-400, 10400)
    plt.ylim(0, 4)


def plot(acceptance, etotal, conf_list, native, nonnative, protein_name, T):

    import numpy as np
    import matplotlib.pyplot as plt

    if protein_name == "mer48A":
        text = "48-mer-A"
    elif protein_name == "mer27":
        text = "27-mer"
    else:
        text = "chignolin"

    steps = len(conf_list)
    stats = int(steps/1000)
    
    # =========================================================================
    # Figure 1 plots the total energy of the system.
    # =========================================================================
    fig1 = plt.figure()
    label1 = "Energy Trajectory of " + text + " at $T$ = " + str(T) + " K"
    plt.scatter(np.arange(0.0, steps, stats), etotal[::stats],
                c='black', s=1, marker='.')
    plt.scatter(np.arange(50 * stats, steps, stats),
                running_mean(etotal[::stats], 51),
                c='blue', s=0.5, marker='.')
    plt.title(label1)
    plt.ylabel("Energy($J/mol$)")
    plt.xlabel("Monte Carlo time (steps)")
    #plt.xlim(-400, 10400)
    plt.ylim(-120, -40)

#    # =========================================================================
#    # Figure 2 plots the contacts map.
#    # =========================================================================
#    fig2 = plt.figure()
#    label2 = "Contacts of " + text + " at $T$ = " + str(T) + " K"
#    plt.scatter(running_mean(native[::stats], 100),
#                running_mean(nonnative[::stats], 100),
#                c='black', s=0.1, marker='.')
#    #plt.xlim(-1,7)
#    #plt.ylim(-1,7)
#    plt.title(label2)
#    plt.ylabel("Non-native Contacts")
#    plt.xlabel("Native Contacts ($Q$)")


    # =========================================================================
    # Figure 4 plots the acceptance ratio trajectory.
    # =========================================================================
    fig4 = plt.figure()
    label4 = "Acceptance Ratio of " + text + " at $T$ = " + str(T) + " K"
    plt.scatter(np.arange(0.0, steps, 1), acceptance,
                c='black', s=1, marker='.')
    plt.scatter(np.arange(50 * stats, steps, stats),
                running_mean(acceptance[::stats], 51),
                c='blue', s=0.5, marker='.')
    plt.title(label4)
    plt.ylabel("% Acceptance")
    plt.xlabel("Monte Carlo time (steps)")
    #plt.xlim(-400, 10400)
    plt.ylim(-0.1, 1.1)

#    # =========================================================================
#    # Figure 5 plots the number of native contacts.
#    # =========================================================================
#    fig5 = plt.figure()
#    label5 = "Native Contacts of " + text + " at $T$ = " + str(T) + " K"
#    plt.scatter(np.arange(0.0, steps, stats), native[::stats],
#                c='black', s=1, marker='.')
#    plt.scatter(np.arange(5000 * stats, steps, stats),
#                running_mean(native[::stats], 5001),
#                c='blue', s=0.0007, marker='.')
#    plt.title(label5)
#    plt.ylabel("$Q$")
#    plt.xlabel("Monte Carlo time (steps)")
#    # plt.xlim(-400, 10400)
#    # plt.ylim(-0.1, 1.1)

    # =========================================================================
    # Figure 6 plots the number of non-native contacts.
    # =========================================================================
    fig6 = plt.figure()
    label6 = "Number of Contacts of " + text + " at $T$ = " + str(T) + " K"
    plt.scatter(np.arange(0.0, steps, stats), nonnative[::stats],
                c='black', s=1, marker='.')
    plt.scatter(np.arange(50 * stats, steps, stats),
                running_mean(nonnative[::stats], 51),
                c='blue', s=0.5, marker='.')
    plt.title(label6)
    plt.ylabel("$Q_n$")
    plt.xlabel("Monte Carlo time (steps)")
    # plt.xlim(-400, 10400)
    # plt.ylim(-0.1, 1.1)

    # =========================================================================
    # Figure 7 plots the standard deviation of the energy trajectory.
    # =========================================================================
    engstd = []
    for i in np.arange(50 * stats, steps-(51*stats), stats):
        # print(i)
        engstd.append(np.std(etotal[(i-50*stats):(i+49*stats):stats]))

    fig7 = plt.figure()
    label7 = "Energy Variance of " + text + " at $T$ = " + str(T) + " K"
    plt.scatter(np.arange(100 * stats, steps-51*stats, stats),
                running_mean(engstd, 51), c='blue', s=0.5, marker='.')
    plt.title(label7)
    plt.ylabel("Energy($J/mol$)")
    plt.xlabel("Monte Carlo time (steps)")
    plt.xlim(-0.1*steps, steps*1.1)
    plt.ylim(0, 8)




#    # =========================================================================
#    # Figure 8 plots the free energy landscape of the lattice protein.
#    # =========================================================================
#    from matplotlib.mlab import griddata
#    import matplotlib.pyplot as plt
#    import numpy as np
#
#    x = np.arange(0,10)
#    y = np.arange(0,25)
#    z = numpy.histogram2d(native, nonnative, bins=[10,25], range=[[0,10],[0,25]])
#    logz = - np.log(z[0] + 1)
#
#    CS = plt.contour(y, x, logz, 15, linewidths=0.5, colors='k')
#    CS = plt.contourf(y, x, logz, 15, vmax=0, vmin=-abs(logz).max())
#    plt.colorbar()  # draw colorbar
#    # plot data points.
#    label8 = "Free Energy Landscape of " + text + " at $T$ = " + T
#    plt.xlim(0,21)
#    plt.title("Free Energy Landscape of 48-mer-A")
#    plt.ylabel("Native Contacts($Q_n$)")
#    plt.xlabel("Non-native Contacts($Q$)")
#    plt.show()

    # =========================================================================
    # Figure 10 plots the Energy distribution.
    # =========================================================================
    fig10 = plt.figure()
    label10 = "Energy distribution of " + text + " at $T$ = " + str(T) + " K"
    bins = np.linspace(-120, -40, 800)
    plt.hist(etotal[::int(stats/100)], bins, log=True)
    plt.title(label10)
    plt.ylabel("Number of counts($N$)")
    plt.xlabel("Energy($J/mol$)")
    #plt.xlim(-400, 10400)
    plt.ylim(1, 5e4)