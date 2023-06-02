"""This file contains a function to calculate the occupation probability of a two-level system
the plot will be animated"""

import numpy as np
from scipy.linalg import expm
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def wavefunction(interaction, energydifference, coefficient1, coefficient2):
    """This function returns the occupation probability of states in a two-level system
    It takes as input:
    :param interaction:the interaction energy between states
    :param energydifference: the energy difference between states
    :param: coefficients 1 and 2: the coefficients specifying the initial population of each states
     it:
        -builds the wavefunction at time 0
        -generates the time-independent hamiltonian matrix
        -defines the basis states phi_1 and phi_2
        -generates wavefunction at all times->list
        -calculates the occupation probability of each state
        -returns the values as float for plotting"""

    # making the time data points
    time = np.linspace(0, 2 * np.pi, 1000)
    # defining the basis states as column vectors
    phi_1 = np.array([[1, 0]]).T
    phi_2 = np.array([[0, 1]]).T
    # Definining the wavefunction at t = 0
    psi_0 = (1 / np.sqrt(coefficient1**2 + coefficient2**2)) * np.array(
        [[coefficient1, coefficient2]]
    ).T
    # Defining the Hamiltonian matrix
    hamiltonian = np.array([[0, interaction], [interaction, energydifference]])
    # generating the wavefunction at all times
    psi_time = [expm(-t * hamiltonian * i).dot(psi_0) for t in time]
    # Defining probability for each state at each time point, and only taking the real part
    probability1 = np.real(
        [(phi_1.T.dot(Psi_t)).dot(Psi_t.conj().T.dot(phi_1)) for Psi_t in psi_time]
    )
    probability2 = np.real(
        [(phi_2.T.dot(Psi_t)).dot(Psi_t.conj().T.dot(phi_2)) for Psi_t in psi_time]
    )
    # The probability will be real but python returns complex where the Im part is zero
    # converting the list of array to a list of floats (needed for plotting)
    prob11 = [float(p1) for p1 in probability1]
    prob22 = [float(p2) for p2 in probability2]
    return prob11, prob22

if __name__ == "__main__":
    i = 1j  # complex i
    INTERACTION = 1  # interaction energy
    ENERGYDIFFERENCE = 0  # energy difference
    timepoints = np.linspace(0, 2 * np.pi, 1000)
    COEFF_1 = 1
    COEFF_2 = 0

    prob1, prob2 = wavefunction(INTERACTION, ENERGYDIFFERENCE, COEFF_1, COEFF_2)

    # creating figure object and subplot to add to the figure
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)

    # making the plot
    line1=ax.plot(timepoints[0], prob1[0],label=f'c1_0={1}, V={1},\u0394={1}')[0]
    line2=ax.plot(timepoints[0], prob2[0],label=f"c2_0={0}")[0]
    ax.set_xlabel("time [a.u.]", fontsize=15)
    ax.set_ylabel("Probability", fontsize=15)
    ax.set(xlim=[0, 7], ylim=[0,1.1])
    ax.legend()
    plt.subplots_adjust(bottom=0.25, top=0.7)

    def update(frames):
        """This function will update the graphs every time the slider is used
        It will make a function call to wavefunction
        It will clear the previous plots and make new plots on the same fig objects
        """
        #unused_value = val  # avoiding pylint error of unused value
        tval = timepoints[:frames]
        yval1 = prob1[:frames]
        yval2 = prob2[:frames]
        # update plots
        line1.set_xdata(tval)
        line1.set_ydata(yval1)
        line2.set_xdata(tval)
        line2.set_ydata(yval2)
        return (line1,line2)



    ani = animation.FuncAnimation(fig=fig, func=update, frames=1000, interval=10)
    #animation.Animation.save(filename="two_levelanim",writer=None,fps=None,dpi=None,codec=None,bitrate=None,extra_args=None,metadata=None,extra_anim=None,savefig_kwargs=None, progress_callback=None)
    ani.save(filename="twolevelanim.gif",writer="Pillow",fps=2)
    plt.show()