def data(direc, saf= True):  
  file = open(direc, 'r')
  wut = file.read()
  forces = wut.split("\n")
  file.close()
  historia_aceleraciones = []
  if saf == True:
    for i, fila in enumerate(forces):
      if i == len(forces)-1: break
      if i == int(SAFFILE["row_duration"]):
        part = fila.split("/")
        duraciones = [abs(float(x)) for x in part[1:]]
      if i == int(SAFFILE["row_accl_max"]):
        part = fila.split("/")
        aceleraciones_maximas = [abs(float(x)) for x in part[1:]]
        a_max = max(aceleraciones_maximas)
        canal = aceleraciones_maximas.index(a_max)
      if i >= int(SAFFILE["row_acl_strt"]):
        part = fila.split(" ")
        aceleraciones = [acel for acel in part if acel.strip()]
        aceleracion_i = float(aceleraciones[canal])
        historia_aceleraciones.append(aceleracion_i)
    time_step = duraciones[canal]/(len(forces) - int(SAFFILE["row_acl_strt"]) - 1)
    return np.array(historia_aceleraciones), time_step
  else:
    for i, aceleracion_i in enumerate(forces):
      historia_aceleraciones.append(float(aceleracion_i))
    return np.array(historia_aceleraciones)

class Show:
  def __init__(self, master):
    self.master = master
    self.frame = tk.Frame(self.master, bg= "white")
    self.frame.place(x= 0, y= 0, relwidth=1, relheight=1)

    self.font_color = "#515151"
    self.title_name = tk.Canvas(self.frame, bd= 0, bg= 'white', highlightthickness= 0, relief= 'flat', height= 600, width= 1000,)
    self.title_name.create_text(315, 70, font=('Open Sans',28), text= "Seismic Response", fill= self.font_color)
    self.title_name.place(x= 0, y= 0 )

    self.open_file = tk.Button(self.frame, text= 'Load Data', bd=0, bg= '#C1F0FF', highlightthickness= 0, relief= 'flat', command= lambda: self.open_data(saf= False),  font=('Open Sans',14), fg= self.font_color)
    self.open_file_saf = tk.Button(self.frame, text= 'Load ASA 2.0 (MEX)', bd=0, bg= '#C1F0FF', highlightthickness= 0, relief= 'flat', command= lambda: self.open_data(saf= True), font=('Open Sans',14), fg= self.font_color)

    self.open_file.place(
      relx= 153.6/1000, rely= 200.64/600,
      relwidth= 258.24/1000, relheight= 63.36/600)
    self.open_file_saf.place(
      relx= 153.6/1000, rely= 349.44/600,
      relwidth= 258.24/1000, relheight= 63.36/600)

  def open_data(self, saf= False):
    self.dir = askopenfilename()
    try:
      self.xi_entry.destroy(); self.xi_label.destroy(); self.mss_label.destroy()
      self.stff_label.destroy(); self.mss_entry.destroy(); self.stff_entry.destroy()
      self.dt_label.destroy(); self.dt_entry.destroy() 
      self.CONTINUAR.destroy(); self.export_.destroy(); self.export_label.destroy()
    except: pass

    try:
      if saf == True:
        self.historia_aceleraciones, self.time_step = data(self.dir, saf= True)
        self.time_step_input = tk.DoubleVar()
        self.time_step_input.set(self.time_step)
      else:
        self.time_step_input = tk.DoubleVar()
        self.historia_aceleraciones = data(self.dir, saf= False)

      self.export  = tk.IntVar()
      self.export_ = tk.Checkbutton(self.frame, bg= "white", variable= self.export, onvalue= 1, offvalue= 0)
      
      self.export_.place(relx=660/1000, rely=340/600)
      
      self.xi_input, self.mss_input, self.stff_input = tk.DoubleVar(), tk.DoubleVar(), tk.DoubleVar()
      
      self.xi_input.set(START['xi'])
      self.mss_input.set(START['m'])
      self.stff_input.set(START['k'])

      self.xi_entry = tk.Entry(self.frame, bd = 2, bg = "#E6E6E6", justify = tk.LEFT, font = ('Consolas', 14), show="", width=20, textvariable= self.xi_input)
      self.mss_entry = tk.Entry(self.frame, bd = 2, bg = "#E6E6E6", justify = tk.LEFT, font = ('Consolas', 14), show="", width=20, textvariable= self.mss_input)
      self.stff_entry = tk.Entry(self.frame, bd = 2, bg = "#E6E6E6", justify = tk.LEFT, font = ('Consolas', 14), show="", width=20, textvariable= self.stff_input)

      self.xi_label = tk.Label(self.frame, text= "ξ (%) =", anchor= "e", fg= self.font_color, bg= "white", font= ("Open Sans", 14), width=10)
      self.mss_label = tk.Label(self.frame, text= "mass (kg) =", anchor= "e", fg= self.font_color, bg= "white", font= ("Open Sans", 14), width=10)
      self.stff_label = tk.Label(self.frame, text= "k (N/m) =", anchor= "e", fg= self.font_color, bg= "white", font= ("Open Sans", 14), width=10)
      self.export_label = tk.Label(self.frame, bg= "white", text= "Export data?", fg= self.font_color, font= ("Open Sans", 14), anchor= "e")

      self.labels_x = 520

      self.xi_label.place(relx=self.labels_x/1000, rely=200/600)
      self.mss_label.place(relx=self.labels_x/1000, rely=235/600)
      self.stff_label.place(relx=self.labels_x/1000, rely=270/600)
      self.export_label.place(relx=self.labels_x/1000, rely=335/600)
      
      self.xi_entry.place(relx=660/1000, rely=200/600)
      self.mss_entry.place(relx=660/1000, rely=235/600)
      self.stff_entry.place(relx=660/1000, rely=270/600)

      self.title_name.delete('baddie')
      if saf == True:
        self.dt_entry = tk.Label(self.frame, bd= 2, text= round(self.time_step,5), bg = "#E6E6E6", relief="sunken", anchor= "w", font = ('Consolas', 14), width=20)
      else:
        self.dt_entry = tk.Entry(self.frame, bd = 2, bg = "#E6E6E6", justify = tk.LEFT, font = ('Consolas', 14), show="", width=20, textvariable= self.time_step_input)
      
      self.dt_label = tk.Label(self.frame, text= "Δt (s) =", anchor= "e", fg= self.font_color, bg= "white", font= ("Open Sans", 14), width=10)
      self.dt_entry.place(relx=660/1000, rely=305/600)
      self.dt_label.place(relx=self.labels_x/1000, rely=305/600)
      
      self.CONTINUAR = tk.Button(self.frame, width= 20, text= "Continue", height= 1, anchor= tk.CENTER, bg="#A6ECA8", font= ("Open Sans", 12), command= self.checker)
      self.CONTINUAR.place(relx=660/1000, rely=400/600)

    except:
      self.title_name.create_text(270, 280, font=('Open Sans',14), text= "Not a valid file. Only 《.txt》", tag= 'baddie', fill="#ff6961")

    return 0

  def checker(self):
    checks = [self.time_step_input, self.mss_input, self.stff_input]
    self.title_name.delete('zeros'); self.title_name.delete('damp')

    for k, p in enumerate(checks):
      if p.get() <= 0:
        messagebox.showwarning("Warning", "No negative numbers or 0")
        return 0

    if self.xi_input.get() > 100 or self.xi_input.get() < 0:
      messagebox.showwarning("Warning", "0< Damp.ratio <100")
      return 0

    seismic_data = {
      "forces": self.historia_aceleraciones * self.mss_input.get() / 100 / 1000,
      "mass": self.mss_input.get(),
      "xi": self.xi_input.get()/100,
      "stiffness": self.stff_input.get(),
      "Dt": self.time_step_input.get(),
      "Export": self.export.get()}

    return Main(seismic_data)

class Main:
  def __init__(self, erthqk_data):
    self.submaster = tk.Tk()
    self.submaster.title("Seismic Response")
    self.submaster.geometry("1280x720")
    self.submaster.minsize(width= 380, height= 380)
    self.submaster.iconbitmap(f'{dir}/assets/pen.ico')
    self.frame = tk.Frame(self.submaster, bg= "white")
    self.submaster.resizable(True,True)
    self.frame.place(x= 0, y= 0, relwidth=1, relheight=1)
    tabs = ttk.Notebook(self.frame)

    self.tab1 = ttk.Frame(tabs)
    self.tab2 = ttk.Frame(tabs)
    self.tab3 = ttk.Frame(tabs)

    tabs.add(self.tab1, text ='History')
    tabs.add(self.tab2, text ='Spectras')
    tabs.add(self.tab3, text ='Mass-Spring')
    tabs.pack(expand = 1, fill ="both")

    self.history(erthqk_data)
    self.spectras(erthqk_data)
    self.animation_earthquake(erthqk_data)

    if erthqk_data["Export"] == 1:
      self.export_data(erthqk_data)

    self.stop_button = tk.Button(self.tab3, bg= "#C1F0FF", text= "Stop", font= ("Open Sans", 14), command= lambda: self.stop_animation())
    self.rsum_button = tk.Button(self.tab3, bg= "#C1F0FF", text= "Play", font= ("Open Sans", 14), command= lambda: self.resume_animation())
    self.stop_button.place(relx= 0.8, rely= 0.5, relwidth= 0.1, relheight= 0.05)
    self.rsum_button.place(relx= 0.8, rely= 0.35, relwidth= 0.1, relheight= 0.05)

  def history(self, erthqk_data):
    """
    History
    """
    h = sr()
    self.r_fig = h.response_fig(erthqk_data, PLOTS)

    self.fig_resp, self.ax_resp = plt.subplots(3, tight_layout= True)
    #self.plot_rspnse_tab1 = tk.Canvas(self.tab2, bd= 0, bg= 'red', highlightthickness= 0, relief= 'flat')
    self.plot_rspnse_tab1 = FigureCanvasTkAgg(self.r_fig, master= self.tab1)
    self.plot_rspnse_tab1.draw()
    self.plot_rspnse_tab1.get_tk_widget().place(x= 0, y= 0, relheight=1, relwidth=1)
    #self.plot_rspnse_tab1.place(x= 0, y= 0, relheight=1, relwidth=0.85)
    if erthqk_data["Export"] == 1:
      u, v, a, t = h.interpol_meth(erthqk_data)
      resp_tabular = {"time" : t,
      "displacement" : u,
      "velocity" : v,
      "acceleration" : a}
      self.resp_data = pd.DataFrame(resp_tabular)
  
  def spectras(self, erthqk_data):
    """
    Spectras
    """
    s = sr()
    solv = s.seismic_spectra(erthqk_data, SPECTRA)
    s_fig = s.seismic_spectra_fig(solv, PLOTS)
    
    self.fig_spec = s_fig
    self.plot_rspnse_tab2 = FigureCanvasTkAgg(self.fig_spec, master= self.tab2)
    self.plot_rspnse_tab2.draw()
    self.plot_rspnse_tab2.get_tk_widget().place(x= 0, y= 0, relheight=1, relwidth=1)

    if erthqk_data["Export"] == 1:
      sp_tabular = {"periods": solv["periods"],
      "maxdisplacement": solv["u_max"],
      "maxvelocities": solv["v_max"],
      "pseudoacceleration": solv["a_max"]}
      self.spc_data = pd.DataFrame(sp_tabular)

  def animation_earthquake(self, seismic_data):
    import matplotlib.pyplot as plt
    import matplotlib.animation as plt_anim
    """
    Animation
    """
    anim_config = ANIMATION
    fig, axes = plt.subplots(1)
    a = sr()
    self.plot_rspnse_tab3 = FigureCanvasTkAgg(fig, master= self.tab3)
    self.plot_rspnse_tab3.get_tk_widget().place(x= 0, y= 0, relheight=1, relwidth=0.7)
    
    u, _, _, time = a.interpol_meth(seismic_data)
    forces = seismic_data['forces']
    Dt = seismic_data['Dt']
    try:
      import seaborn
      seaborn.despine(fig, top= True, right= True, left= True)
    except: pass
    axes.set_xlim(min(u),max(u)); axes.set_ylim(0,2.5)
    axes.set_xlabel("u [dm]")
    axes.get_yaxis().set_visible(False)

    mass_dot_size = int(anim_config['mass_dot_size'])
    mass_dot_color = anim_config['mass_dot_color']

    def animate(i):
      height = np.linspace(0,2,1000)
      if u[i] != 0:
        seisEI = forces[i] / u[i] * (2 * 2 ** 3 - 3 * 2 ** 2 * 0 + 0 ** 3)
        if seisEI != 0:
          colux = forces[i] / seisEI * (2 * 2 ** 3 - 3 * 2 ** 2 * np.flip(height) + np.flip(height) ** 3)
        else: 
          colux = np.zeros(height.shape)
      else: 
        colux = np.zeros(height.shape)
      
      line2, = axes.plot(colux/10, height, color= "black", lw= 2)
      line, = axes.plot([u[i]/10], [2], marker= "o", markersize= mass_dot_size, color= mass_dot_color)
      
      tm = axes.text(0, 2.4, 'Time = %.1fs' % time[i])
      return line2, line, tm
      
    self.ani = plt_anim.FuncAnimation(fig, animate, frames= len(u), interval= Dt, blit= True, repeat= True)
    self.plot_rspnse_tab3.draw()

  def stop_animation(self):
    return self.ani.event_source.stop()

  def resume_animation(self):
    return self.ani.event_source.start()

  def export_data(self,erthqk_data):
    name = f'{erthqk_data["mass"]}-{erthqk_data["stiffness"]}-{erthqk_data["xi"]}-{erthqk_data["Dt"]}'
    self.spc_data.to_csv(f'{dir}\{name}_spectra.txt', sep= '\t')
    self.resp_data.to_csv(f'{dir}\{name}_response.txt', sep= '\t')
    return 0

if __name__ == '__main__':
  from tkinter.filedialog import askopenfilename
  from configparser import ConfigParser
  from tkinter import ttk, messagebox
  from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
  import matplotlib.pyplot as plt
  import tkinter as tk, numpy as np, pandas as pd, os
  
  dir = os.path.dirname(os.path.realpath(__file__))
  config = ConfigParser()
  config.read(f"{dir}/config.ini")
  START = config['START']
  SAFFILE = config['SAFFILE']
  PLOTS = config['PLOTS']
  ANIMATION = config['ANIMATION']
  SPECTRA = config['SPECTRA']

  from seismic_response import Seismic as sr

  root = tk.Tk()
  root.title("Seismic Response")
  root.geometry("1000x600")
  root.iconbitmap(f'{dir}/assets/pen.ico')
  root.resizable(False, False)

  a = Show(root)

  root.mainloop()