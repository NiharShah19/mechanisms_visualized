# import the necessary modules
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import TextBox, CheckButtons

# create an initial figure, add an axes and add lines into this axes
fig,ax = plt.subplots()
plt.subplots_adjust(left=0.3,bottom=0.3)
ax.set_xlabel('Velocity(m/s)')
ax.set_ylabel('Thrust Required(N)')
x = np.array([0])
line1, = ax.plot(x, x, color='blue')
line2, = ax.plot(x, x, color='red')
line3, = ax.plot(x, x, color='black')
lines = [line1, line2, line3]

# create a text box widget
ax1 = fig.add_axes([0.4,0.1,0.2,0.1])
rho = TextBox(ax=ax1,label='density(kg/$m^3$)',hovercolor='yellow',initial='1.15')
density = '1.15'

# define the sea level density
rho_sl = 1.225

# create a class named 'jetAircrafts', define the attributes and methods
# which will be available to all the instances


class jet_aircrafts():

    # to set the initial values of attributes
    def __init__(self, b, s, w, T_sl_max, C_do, e):
        self.wing_span = b
        self.wing_area = s
        self.weight = w
        self.max_thrust = T_sl_max*2
        self.C_do = C_do
        self.e = e
        self.k = np.pi*self.e*(self.wing_span**2/self.wing_area)
        self.vel_kmph = None
        self.p_4 = None
        self.p_2 = None
        self.p_0 = None
        self.p = None
        self.solnT = None
        self.solnA = None
        self.minVel = None
        self.maxVel = None

    # converting result from mps to kmph
    def mps_to_kmph(self, vel_mps):
        self.vel_kmph = (10.0/36.0)*vel_mps
        return self.vel_kmph

    # function to obtain the roots by solving the polynomial
    def get_roots(self, rho_al):
        self.p_4 = 0.25 * (rho_al ** 2) * self.k * self.C_do * (self.wing_area ** 2)
        self.p_2 = -0.5 * rho_al * self.wing_area * self.k * ((rho_al / rho_sl) * self.max_thrust)
        self.p_0 = self.weight ** 2

        self.p = np.array([self.p_4, 0, self.p_2, 0, self.p_0])
        self.solnT = np.roots(self.p)
        self.solnA = self.solnT[self.solnT[np.isreal(self.solnT)] > 0]

    # function to find minimum and maximum velocity(positive)
    def min_max(self):

        self.minVel = self.solnA[0]
        for sol in self.solnA:
            if sol < self.minVel:
                self.minVel = sol
            else:
                self.minVel = self.minVel

        self.maxVel = self.solnA[0]
        for sol in self.solnA:
            if sol > self.maxVel:
                self.maxVel = sol
            else:
                self.maxVel = self.maxVel


# instantiation
cessna = jet_aircrafts(16.25, 29.54, 88176.75, 16242.5, 0.02, 0.81)
boeing_737max = jet_aircrafts(35.92, 127, 80000*9.8, 130000, 0.019, 0.81)
airbus_a380 = jet_aircrafts(79.75, 845, 575000*9.8, 348000*2, 0.0265, 0.81)

# additional text boxes to print the maximum and minimum velocity in case only one
# object is instantiated
        # ax2 = fig.add_axes([0.5,0.1,0.1,0.1])
        # minV = TextBox(ax=ax2,label='min_vel',hovercolor='yellow')
        #
        # ax3 = fig.add_axes([0.8,0.1,0.1,0.1])
        # maxV = TextBox(ax=ax3,label='max_vel',hovercolor='yellow')

# create a checkButtons widget
ax4 = fig.add_axes([0.02, 0.4, 0.2, 0.3])
label = ['cessna_citationIII', 'boeing_737', 'airbus_a380']
chkbx = CheckButtons(ax4, label)

#######################################################
# defining some lists and dictionaries to be used later
my_dict = {'cessna':0, 'boeing':1, 'airbus':2}
objs = [cessna, boeing_737max, airbus_a380]
max_xlim = []
max_ylim = []
# objs_selected = []
all_indices = [0, 1, 2]
#######################################################

# function to redraw the plot based on the data obtained from 'get_roots' and 'min_max' function


def my_fun(rho_al,new_ind):
    for i in all_indices:

        if i in new_ind:

            for ind in new_ind:
                global max_xlim, max_ylim
                max_xlim = []
                max_ylim = []

                objs[ind].get_roots(float(rho_al))
                objs[ind].min_max()
                v = np.linspace(objs[ind].minVel,objs[ind].maxVel,100)
                q_inf = 0.5*float(rho_al)*v**2
                T_req = (q_inf*objs[ind].wing_area*objs[ind].C_do) + \
                        (objs[ind].weight**2/(q_inf*objs[ind].wing_area*objs[ind].k))

                # minV.set_val(np.round(objs[ind].minVel,2))
                # maxV.set_val(np.round(objs[ind].maxVel,2))
                lines[ind].set_xdata(v)
                lines[ind].set_ydata(T_req)
                max_xlim.append(objs[ind].maxVel)
                max_ylim.append(T_req[99])
                # set the label property of line when the corresponding checkbox is selected
                lines[ind].set_label(list(my_dict.keys())[list(my_dict.values()).index(ind)])

        else:
            lines[i].set_xdata(x)
            lines[i].set_ydata(x)
            lines[i].set_label('')


    ax.set_xlim(0, np.max(np.array(max_xlim)) + 10)
    ax.set_ylim(0, np.max(np.array(max_ylim))+10000)
    ax.legend()
    plt.draw()

# function to call the 'my_func' function when a change is made to the 'rho' textbox object


def select_density(den):
    global density
    density = den
    my_tup = chkbx.get_status()
    true_ind = [i for i, x in enumerate(my_tup) if x]
    my_fun(density, true_ind)

    # print(true_ind)

# function to call the 'my_func' function when a change is made to the 'chkbx' checkButton object


def select_aircraft(val):
    # i = my_dict[val]
    # objs_selected.append(i)
    my_tup = chkbx.get_status()
    true_ind = [i for i, x in enumerate(my_tup) if x]
    my_fun(density, true_ind)
    

# call to the function 'select_density' when a change is made in density


rho.on_submit(select_density)

# call to the function 'select_aircraft' when some checkButton is selected


chkbx.on_clicked(select_aircraft)
plt.show()








