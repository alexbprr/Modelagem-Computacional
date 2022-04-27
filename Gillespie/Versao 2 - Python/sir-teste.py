import gillespie
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

if __name__ == '__main__':
    N = 500  # whole population
    beta = 3  # transmission rate
    gamma = 0.3  # recovery rate
    t = 15  # duration

    I0 = 3
    initials = [N - I0, I0, 0]  # S, I, R

    propensities = [lambda s, i, r: (beta * s * i) / N,  # S -> I, Propensity: b * S(t) * I(t) / N
                    lambda s, i, r: gamma * i]  # I -> R Propensity: g * I(t)

    stoichiometry = [[-1, 1, 0],  # S -> I, Population change: S-1, I+1, R+0
                     [0, -1, 1]]  # I -> R Population change: S+0, I-1, R+1

    pdf = PdfPages('SIR.pdf')
    pdf2 = PdfPages('SIR-all.pdf')

    figAll = plt.figure(figsize=(6,4)) 
    plt.title("SIR epidemic model")
    plt.xlabel("Days")
    plt.ylabel("Population")    
    ax = figAll.gca()

    nRuns = 20 
    k = 0
    while k < nRuns:                  
        t = 15
        t, SIR = gillespie.simulate(initials, propensities, stoichiometry, t)
        S, I, R = zip(*SIR)

        fig = plt.figure(figsize=(6,4))
        plt.plot(t, S, label="susceptible")
        plt.plot(t, I, label="infected")
        plt.plot(t, R, label="recovered")

        plt.title("SIR epidemic model")
        plt.xlabel("Days")
        plt.ylabel("Population")
        plt.legend()
        #plt.savefig("plots/SIR" + str(k) + ".png")
        pdf.savefig(fig)

        ax.plot(t, S, color="blue", label="susceptible")
        ax.plot(t, I, color="orange", label="infected")
        ax.plot(t, R, color="green", label="recovered")

        k += 1

    pdf2.savefig(figAll)
    pdf2.close()
    pdf.close()
