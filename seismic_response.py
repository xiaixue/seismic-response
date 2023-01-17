import numpy as np
from math import sin, cos, e

def data(dir):  
  file = open(dir, 'r')
  wut = file.read()
  forces = wut.split("\n")
  file.close()
  historia_aceleraciones = []
  for i, fila in enumerate(forces):
    if i == len(forces)-1: break
    if i == 69:
      part = fila.split("/")
      duraciones = [abs(float(x)) for x in part[1:]]
    if i == 73:
      part = fila.split("/")
      aceleraciones_maximas = [abs(float(x)) for x in part[1:]]
      a_max = max(aceleraciones_maximas)
      canal = aceleraciones_maximas.index(a_max)
    if i >= 109:
      part = fila.split(" ")
      aceleraciones = [acel for acel in part if acel.strip()]
      aceleracion_i = float(aceleraciones[canal])
      historia_aceleraciones.append(aceleracion_i)
  time_step = duraciones[canal]/(len(forces) - 111)
  return historia_aceleraciones, time_step

class Seismic:
  def interpol_meth(self, seismic_data):
    m = seismic_data["mass"]
    k = seismic_data["stiffness"]
    xi = seismic_data["xi"]
    forces = seismic_data["forces"]
    Dt = seismic_data["Dt"]
    v0 = 0 # forces[0]*1000 * Dt / m
    u0 = 0 # forces[0]*1000 / k
    omega = (k / m) ** 0.5
    omega_d = omega * (1 - xi ** 2) ** 0.5
    u, u_ = [], []

    A_ = (e ** (-xi * omega * Dt)) * ( ( xi/((1-xi**2)** 0.5)) * sin(omega_d*Dt) + cos(omega_d*Dt))
    B_ = (e ** (-xi * omega * Dt)) * ( (1/omega_d) * sin(omega_d*Dt))
    C_ = (k**(-1)) * (2*xi/omega/Dt + (e ** (-xi*omega*Dt)) * ( ( (1 - 2 * xi * xi)/omega/Dt - xi/ (( 1 - xi ) ** 0.5) ) * sin(omega*Dt)- (1 + 2*xi/omega/Dt) * cos(omega_d*Dt))) * 10e4
    D_ = (k ** -1) * (1 - 2* xi /omega/Dt + (e **(-xi*omega*Dt)) * ( (2*xi*xi -1) * sin(omega_d*Dt)/omega/Dt + 2*xi*cos(omega_d*Dt)/omega/Dt)) * 10e4
    A__ = (-e ** (-xi*omega*Dt)) * (omega*sin(omega_d*Dt)/((1-xi*xi)**0.5))
    B__ = (e ** (-xi*omega*Dt)) * (cos(omega_d*Dt) - xi * sin(omega_d*Dt) / ((1-xi*xi)**0.5))
    C__ = (k**-1) * (-1/Dt + e ** (-xi*omega*Dt) * ((omega_d/((1-xi*xi)**0.5) + xi / Dt / ((1-xi*xi)**0.5)) * sin(omega_d*Dt) + cos(omega_d*Dt)/Dt)) * 10e4
    D__ = 1/k/Dt * (1 - e ** (-xi * omega * Dt) * (xi * sin(omega_d*Dt)/((1-xi*xi)**0.5) + cos(omega_d*Dt))) * 10e4
    u_.append(u0); u.append(v0)
    
    for i in range(len(forces)):
      if i == 0: continue
      else:
        velocity = A__ * u0 + B__ * v0 + C__ * forces[i-1] + D__ * forces[i]
        motion = A_ * u0 + B_ * v0 + C_ * forces[i-1] + D_ * forces[i]
      u0 = motion; v0 = velocity
      u_.append(velocity); u.append(motion)
    
    u__ = []
    for l in range(len(u)):
      u__i = - omega ** 2 * u[l] - 2 * xi * omega * u_[l] + forces[l] * 1000 / m / 100
      u__.append(u__i)

    return u, u_, u__, np.linspace(0.0, len(forces) * Dt, len(forces))

  def seismic_spectra(self, seismic_data, spectras_config):
    T_max, points = float(spectras_config["t_max"]), int(spectras_config["spctrum_pts"])
    self.xi = seismic_data["xi"]

    periodos = np.linspace(0.01, T_max, points)
    u_maxes, u_1_maxes, u_2_maxes = [],[],[]

    for T in periodos:
      stiffness = ( 2 * np.pi / T ) ** 2 * seismic_data["mass"] 
      seismic_data["stiffness"] = stiffness
      u, u_, u__, _ = self.interpol_meth(seismic_data)

      if max(u__) > abs(min(u__)): u_2_max = max(u__)
      else: u_2_max = abs(min(u__))
      u_2_maxes.append(u_2_max)

      if max(u_) > abs(min(u_)): u_1_max = max(u_)
      else: u_1_max = abs(min(u_))
      u_1_maxes.append(u_1_max)

      if max(u) > abs(min(u)): u_max = max(u)
      else: u_max = abs(min(u))
      u_maxes.append(u_max)

    return {"periods": periodos, "u_max": u_maxes, "v_max": u_1_maxes, "a_max": u_2_maxes}
  
  def seismic_spectra_fig(self, spectra_data, plot_config):
    import matplotlib.pyplot as plt
    periodos = spectra_data["periods"]
    u_maxes = spectra_data["u_max"]
    u_1_maxes = spectra_data["v_max"]
    u_2_maxes = spectra_data["a_max"]

    gridstyle = plot_config["gridstyle"]
    lwd = int(plot_config["lw"])
    hstry_a = plot_config["hstry_a"]
    hstry_d = plot_config["hstry_d"]
    hstry_v = plot_config["hstry_v"]
    fig_spec, ax = plt.subplots(ncols= 3, tight_layout=True)

    ax[0].plot(periodos, np.array(u_maxes), color= hstry_d, label= f'ξ(%)= {self.xi*100} %', lw = lwd)
    ax[1].plot(periodos, np.array(u_1_maxes), color= hstry_v, label= f'ξ(%)= {self.xi*100} %', lw = lwd)
    ax[2].plot(periodos, np.array(u_2_maxes), color= hstry_a, label= f'ξ(%)= {self.xi*100} %', lw = lwd)

    ax[0].set_xlim(0, max(periodos)); ax[1].set_xlim(0, max(periodos)); ax[2].set_xlim(0, max(periodos))
    ax[0].set_ylim(0, max(u_maxes)*1.1); ax[1].set_ylim(0, max(u_1_maxes)*1.1); ax[2].set_ylim(0, max(u_2_maxes)*1.1)
    ax[0].legend(); ax[1].legend(); ax[2].legend()
    ax[0].set_xlabel("T (s)"); ax[1].set_xlabel("T (s)"); ax[2].set_xlabel("T (s)")
    ax[0].set_ylabel(r"$u_{max}$ [cm]"); ax[1].set_ylabel(r"$\dot{u}_{max}$ [cm/s]"); ax[2].set_ylabel(r"pseudo $\ddot{u}_{max}$ [gal]")
    ax[0].grid(linestyle=gridstyle); ax[1].grid(linestyle=gridstyle); ax[2].grid(linestyle=gridstyle)  

    return fig_spec
  
  def response_fig(self, seismic_data, plot_config):
    import matplotlib.pyplot as plt
    disp, velo, acel, time = self.interpol_meth(seismic_data)
    time_step = seismic_data["Dt"]

    gridstyle = plot_config["gridstyle"]
    lwd = int(plot_config["lw"])
    hstry_a = plot_config["hstry_a"]
    hstry_d = plot_config["hstry_d"]
    hstry_v = plot_config["hstry_v"]

    fig_resp, ax = plt.subplots(nrows= 3, tight_layout=True)
    
    ax[0].plot(time , np.array(disp), color= hstry_d, linewidth=lwd)
    ax[1].plot(time , np.array(velo), color= hstry_v, linewidth=lwd)
    ax[2].plot(time , np.array(acel), color= hstry_a, linewidth=lwd)

    ax[0].set_xlim(0,len(disp) * time_step - time_step)
    ax[1].set_xlim(0,len(velo) * time_step - time_step)
    ax[2].set_xlim(0,len(acel) * time_step - time_step)

    ax[0].set_ylim(min(disp) * 1.1, max(disp) * 1.1)
    ax[1].set_ylim(min(velo) * 1.1, max(velo) * 1.1)
    ax[2].set_ylim(min(acel) * 1.1, max(acel) * 1.1)

    ax[0].set_ylabel("$u(t)$ [cm]"); ax[1].set_ylabel("$\dot{u}(t)$ [cm/s]"); ax[2].set_ylabel("$\ddot{u}(t)$ [gal]")
    ax[2].set_xlabel("$t$ [s]")

    ax[0].grid(linestyle=gridstyle); ax[1].grid(linestyle=gridstyle); ax[2].grid(linestyle=gridstyle)

    return fig_resp

  def move_backup(self, seismic_data):
    import matplotlib.pyplot as plt
    import matplotlib.animation as anim
    u, _, _, _ = self.interpol_meth(seismic_data)
    forces = seismic_data['forces']

    time = np.arange(0.0,len(u)*time_step, time_step)
    fig, axes = plt.subplots(1)
    
    for i, j in enumerate(time):
      height = np.linspace(0,2,1000)
      
      if u[i] != 0:
        seisEI = forces[i] / u[i] * (2 * 2 ** 3 - 3 * 2 ** 2 * 0 + 0 ** 3)
        if seisEI != 0:
          colux = forces[i] / seisEI * (2 * 2 ** 3 - 3 * 2 ** 2 * np.flip(height) + np.flip(height) ** 3)
        else: 
          colux = np.zeros(height.shape)
      else: 
        colux = np.zeros(height.shape)

      axes.plot(colux/10, height, color= "black",lw= 2)
      axes.plot([u[i]/10], [2], marker= "o", markersize= mass_dot_size, color=mass_dot_color)

      axes.set_xlim(min(u),max(u))
      axes.set_ylim(0,2.5)
      axes.set_xlabel("Desplazamiento [dm]")
      tm = axes.text(0, 2.4, 'Time = %.1fs' % j)
      axes.get_yaxis().set_visible(False)
      plt.pause(0.000000000001)
      axes.cla()
    return

if '__main__' == __name__:
  dit_1 = './Seismic Response/Sismos historicos/Oaxaca 73/OAXM7308_199.txt'
  dir = dit_1
  
  m = 45594   # kg
  k = 1800000 # N / m
  xi = 0.05

  gridstyle = '--'
  lw = 1
  color_spectra_a = '#000'
  hstry_a = '#215968'
  hstry_d = '#93CDDD'
  hstry_v = '#A7C7E7'
  mass_dot_size = 60
  mass_dot_color = '#93CDDD'

  historia_aceleraciones, time_step = data(dir)

  seismic_data = {
    'mass': m,
    'stiffness': k,
    'xi': xi,
    'Dt': time_step,
    'forces': historia_aceleraciones}

  a = Seismic()
  forces = np.array(historia_aceleraciones) * m / 100 / 1000
  sol, sol1, u__, t = a.interpol_meth(seismic_data) 