import numpy as np
import matplotlib.pyplot as plt

def probability_partner_inside(pmax, alpha, x0, mfp):
    d = (x0 + np.cos(alpha)*mfp)/np.sin(alpha)
    if d > pmax:
        if x0 + mfp*np.cos(alpha) - pmax*np.sin(alpha) > 0.:
            return 1.
        else:
            return 0.
    A = pmax*(pmax*np.arccos(d/pmax) - d*np.sqrt(1. - d**2/pmax**2))
    if np.isnan(A):
        return 0
    A_max = np.pi*pmax**2
    return (A_max - A)/A_max

def integral(R, x):
    return 0.5*(x*np.sqrt(R**2 - x**2) + R**2*np.arctan2(x, np.sqrt(R**2 - x**2)))

def depth_distribution(pmax, alpha, x0, mfp, x):
    if x < 0.: return 0.
    d = (x0 + np.cos(alpha)*mfp)/np.sin(alpha)
    f = 2.*np.sqrt(pmax**2 - (x - d)**2)*np.cos(alpha)

    if np.isnan(f):
        return 0.
    else:
        return f

def main():
    mfp = 1.
    pmax = mfp/np.sqrt(2.*np.pi)
    alpha_deg = 60
    alpha = alpha_deg*np.pi/180.
    x0 = np.linspace(-2, 0., 10000)

    example(pmax, alpha, -2.*pmax, mfp)

    plt.figure(1)
    p_min = [[], [], [], []]
    p_max = [[], [], [], []]

    for k in (0, 1, 2, 3):
        for x in x0:
            p_min[k].append(probability_partner_inside(pmax + k*pmax, alpha, x, 0.))
            p_max[k].append(probability_partner_inside(pmax + k*pmax, alpha, x, mfp))


    for p in p_min, p_max:
        p[1] = 1 - (1 - np.array(p[0]))*(1 - np.array(p[1]))
        p[2] = 1 - (1 - np.array(p[0]))*(1 - np.array(p[1]))*(1 - np.array(p[2]))
        p[3] = 1 - (1 - np.array(p[0]))*(1 - np.array(p[1]))*(1 - np.array(p[2]))*(1 - np.array(p[3]))

    plots = []
    for k in (0, 1, 2, 3):
        p_min_plot = plt.plot(x0, p_min[k], linestyle='--')
        plt.plot(x0, p_max[k], color=p_min_plot[0].get_color(), linestyle='--')
        fill_plot = plt.fill_between(x0, p_min[k], p_max[k], color=p_min_plot[0].get_color(), alpha = 0.25)
        plots.append(fill_plot)

    tridyn_line = plt.plot(-2.*pmax*np.ones(10), np.linspace(0., 1., 10), linestyle='--')
    plots.append(tridyn_line[0])

    proposal_line = plt.plot(-mfp*np.cos(alpha)*np.ones(10), np.linspace(0., 1., 10), linestyle='--')
    plots.append(proposal_line[0])

    plt.legend(plots, ('0 weak collisions', '1', '2', '3', 'TRIDYN starting depth', 'x0 = -mfp*cos(alpha)'))

    plt.title(f'BCA Probability of finding first collision partner for angle={alpha_deg}')
    plt.xlabel('initial depth x0 [mfp]')
    plt.ylabel('P(x)')

    plt.show()

def example(pmax, alpha, x0, mfp):
    y0 = 0.

    ca = np.cos(alpha)
    cb = np.sin(alpha)
    cg = 0.

    x1 = x0 + ca*mfp
    y1 = cb*mfp

    #let x_recoil: f64 = x + mfp*ca - impact_parameter*cphi*sa;
    #let y_recoil: f64 = y + mfp*cb - impact_parameter*(sphi*cg - cphi*cb*ca)/sa;
    cphi = 1.

    xa = x0 + mfp*ca - pmax*cphi*np.sin(alpha)
    ya = y0 + mfp*cb - pmax*(-cphi*cb*ca)/np.sin(alpha)
    xb = x0 + mfp*ca + pmax*cphi*np.sin(alpha)
    yb = y0 + mfp*cb - pmax*(cphi*cb*ca)/np.sin(alpha)

    p_zero = x1/np.sin(alpha)
    xz = 0.
    yz = y0 + mfp*cb + p_zero*(cphi*cb*ca)/np.sin(alpha)

    d_zero = (x0 + np.cos(alpha)*mfp)/np.sin(alpha)
    h = d_zero - pmax
    c = 2*pmax*np.sqrt(1 - (d_zero/pmax)**2)
    s = np.arcsin(c/(h + c**2/4./h))*(h + c**2/4./h)
    A = pmax*(pmax*np.arccos(d_zero/pmax) - d_zero*np.sqrt(1. - d_zero**2/pmax**2))
    A_max = np.pi*pmax**2

    print((A_max - A)/(A_max))

    plt.scatter(y0, x0)
    plt.scatter(y1, x1)
    plt.plot([y0, y1], [x0, x1])
    plt.plot([ya, yb], [xa, xb])
    plt.scatter(yz, xz)
    plt.plot(np.linspace(-2, 2, 10), np.zeros(10))
    plt.axis('square')
    plt.show()

if __name__ == '__main__':
    main()
