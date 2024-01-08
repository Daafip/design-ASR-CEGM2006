# import the necessary packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
plt.rcParams['figure.figsize'] = (5, 3) # set default figure size
import flopy as fp  # import flopy and call it fp
import glob
import imageio
import IPython
from flopy.mf6.mfbase import VerbosityLevel


k_lst = [10, 20, 30, 40]
npor_lst = [0.2, 0.3, 0.4, 0.5]
n_run = 0
for k in k_lst:
    for npor in npor_lst:
        print(n_run,end="\n")
        Q_d = 40_000
        Q_tot = Q_d * 1.25 ## Change this later as we actually want to produce this, currentlty we still have 80% loss
        d_extrating = 62
        d_injecting = 365 - d_extrating ## change later for days with water excess
        # print(f'Need to pump {Q_tot/d_extrating:.2f}m^3/d to full fill demand') 
        
        # domain size and boundary conditions
        R = 200 # length of domain, m
        hR = 0 # head at r=R
        
        # aquifer parameters
        # k = 40 # hydraulic conductivity, m/d -> ############ 10! 
        H = 20 # aquifer thickness, m - fixed
        # npor = 0.35 # porosity, generally 0.25 - 0.5  ############# 0.25
        
        # flow
        
        Q = Q_tot/d_extrating # extraction rate, m^3/d 
        Q_in = Q_tot / (d_injecting)# injection rate m^3/d
        
        # transport
        alphaL = 0.5 # longitudinal dispersivity in horizontal direction, m - # something to check!! -> slides guest lecture
        alphaT = alphaL / 10 # transverse dispersivity in horizontal direction, m
        diffusion_coef = 0 # diffusion coeffcient
        
        # concentration
        cs = 30 # initial concentration, kg/m^3 (=g/L)
        cf = 0 # concentration injected water, kg/m^3 (=g/L)
        
        # space discretization
        delr = 0.2 # length of cell along row (in x-direction), m
        delc = 1 # width of cells normal to plane of flow (in y-direction), m
        z = [0, -H] # top and bottom(s) of layers
        nlay = 1 # number of layers
        nrow = 1 # number of rows
        ncol = round(R / delr) # number of columns
        rw = 0.2 # radius of well, m
        r = rw + np.cumsum(delr * np.ones(ncol)) - 0.5 * delr # radial coordinates of centers of cells, m
        
        # convert parameters for radial flow following procedure of Langevin, 2008
        theta = 2 * np.pi # angle for full circle
        krad = k * r * theta # convert k for radial flow
        nporrad = npor * r * theta # convert porosity for radial flow
        
        # time discretization
        ########## injection
        tin = d_injecting # injection time, d - rest of the year - maybe change later
        delt = 0.1 # time step, d
        nstepin = round(tin / delt) # computed number of steps during injection, integer
        
        ######### extraction
        tout = d_extrating # extraction time, d
        delt = 0.1 # time step, d
        nstepout = round(tout / delt) # computed number of steps during extraction, integer
        
        # model name and workspace
        modelname = f'modelrad_{k}_{str(npor).split(".")[-1]}' # name of model
        gwfname = modelname + 'f' # name of flow model
        gwtname = modelname + 't' # name of transport model
        modelws = './' + modelname # model workspace to be used
        
        # simulation
        sim = fp.mf6.MFSimulation(sim_name=modelname, # name of simulation
                                  version='mf6', # version of MODFLOW
                                  exe_name='../bin/mf6', # absolute path to MODFLOW executable
                                  sim_ws=modelws, # path to workspace where all files are stored
                                 )
        
        # time discretization
        tdis = fp.mf6.ModflowTdis(simulation=sim, # add to the simulation called sim (defined above)
                                  time_units="DAYS", 
                                  nper=2, # number of stress periods 
                                  perioddata=[[tin, nstepin, 1],
                                              [tout, nstepout, 1]], # period length, number of steps, timestep multiplier
                                 )
        
        # groundwater flow model
        gwf = fp.mf6.ModflowGwf(simulation=sim, # add to simulation called sim
                                modelname=gwfname, # name of gwf model
                                save_flows=True, # make sure all flows are stored in binary output file
                               )
        
        # iterative model solver
        gwf_ims  = fp.mf6.ModflowIms(simulation=sim, # add to simulation called sim
                                     filename=gwf.name + '.ims', # file name to store ims
                                     linear_acceleration="BICGSTAB", # use BIConjuGantGradientSTABalized method
                                    )                                                                                                
        # register solver
        sim.register_ims_package(solution_file=gwf_ims, # name of iterative model solver instance
                                 model_list=[gwf.name], # list with name of groundwater flow model
                                )   
        
        # discretization
        gwf_dis = fp.mf6.ModflowGwfdis(model=gwf, # add to groundwater flow model called gwf
                                       nlay=nlay, 
                                       nrow=nrow, 
                                       ncol=ncol, 
                                       delr=delr, 
                                       delc=delc, 
                                       top=z[0], 
                                       botm=z[1:], 
                                      )
        
        # aquifer properties
        gwf_npf  = fp.mf6.ModflowGwfnpf(model=gwf, 
                                        k=krad, # horizontal k value
                                        alternative_cell_averaging="LOGARITHMIC", # logarithmic averaging
                                        save_flows=True, # save the flow for all cells
                                       )
            
        # initial condition
        gwf_ic = fp.mf6.ModflowGwfic(model=gwf, 
                                     strt=hR, # initial head used for iterative solution
                                    )
        
        # wells  ############## change when we inject
        wellin = [[(0, 0, 0),  Q_in, cf]]   # [(layer, row, col), U, concentration]
        wellout = [[(0, 0, 0),  -Q, cf]] # specified concentration is not used, but must be specified 
        wel_spd = {0: wellin, 1: wellout} # stress period data for periods 0 and 1
        gwf_wel = fp.mf6.ModflowGwfwel(model=gwf, 
                                       stress_period_data=wel_spd, 
                                       auxiliary=['CONCENTRATION'],
                                       pname='WEL1', # package name
                                      )
        
        # constant head 
        chd0 = [[(0,  0,  ncol-1), hR, cs]] # [(layer, row, col), head, concentration]
        chd_spd  = {0: chd0}    # Stress period data
        gwf_chd = fp.mf6.ModflowGwfchd(model=gwf, 
                                       stress_period_data=chd_spd, 
                                       auxiliary=['CONCENTRATION'],
                                       pname='CHD1', # package name
                                      )
            
        # output control
        oc = fp.mf6.ModflowGwfoc(model=gwf, 
                                 saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")], # what to save
                                 budget_filerecord=f"{gwfname}.cbc", # file name where all budget output is stored
                                 head_filerecord=f"{gwfname}.hds", # file name where all head output is stored
                                )
        
        # groundwater transport model
        gwt = fp.mf6.ModflowGwt(simulation=sim, 
                                modelname=gwtname, # name of groundwater transport model
                               )
        
        # iterative model solver
        gwt_ims  = fp.mf6.ModflowIms(simulation=sim,
                                     filename=gwt.name + '.ims', # must be different than file name of gwf model ims
                                     linear_acceleration="BICGSTAB",
                                    ) 
        sim.register_ims_package(solution_file=gwt_ims, 
                                 model_list=[gwt.name],
                                )
        
        # discretization
        gwt_dis = fp.mf6.ModflowGwtdis(model=gwt, 
                                       nlay=nlay, 
                                       nrow=nrow, 
                                       ncol=ncol, 
                                       delr=delr, 
                                       delc=delc, 
                                       top=z[0], 
                                       botm=z[1:], 
                                      )
        
        # mobile storage and transfer
        gwt_sto = fp.mf6.ModflowGwtmst(model=gwt, 
                                       porosity=nporrad, # porosity
                                       save_flows=True,
                                      )
        
        # initial condition
        gwt_ic = fp.mf6.ModflowGwtic(model=gwt, 
                                     strt=cs, # initial concentration
                                    ) 
        
        # source sink mixing
        sourcelist = [("WEL1", "AUX", "CONCENTRATION"), ("CHD1", "AUX", "CONCENTRATION")]
        ssm = fp.mf6.ModflowGwtssm(model=gwt, 
                                   sources=sourcelist, 
                                   save_flows=True,
                                   pname='SSM1', 
                                  )
        
        # advection
        adv = fp.mf6.ModflowGwtadv(model=gwt,  
                                   scheme="TVD",  # use Total Variation Diminishing (TVD)
                                   pname='ADV1',
                                  )
        
        # dispersion
        dsp = fp.mf6.ModflowGwtdsp(model=gwt, 
                                   alh=alphaL,
                                   ath1=alphaT, 
                                   diffc=diffusion_coef,
                                   pname='DSP1', 
                                  )
        
        # output control
        oc = fp.mf6.ModflowGwtoc(model=gwt,
                                 saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")], # what to save
                                 budget_filerecord=f"{gwtname}.cbc", # file name where all budget output is stored
                                 concentration_filerecord=f"{gwtname}.ucn", # file name where all concentration output is stored
                                )
        
        fp.mf6.ModflowGwfgwt(simulation=sim, 
                             exgtype="GWF6-GWT6", 
                             exgmnamea=gwf.name , 
                             exgmnameb=gwt.name , 
                             filename=f"{modelname}.gwfgwt",
                            );
        
        sim.write_simulation(silent=True)
        success, _ = sim.run_simulation(silent=True) 
        if success == 1:
            print('Model solved successfully',end=" ")
        else:
            print('Solve failed')
        
        cobj = gwt.output.concentration() # get handle to binary concentration file
        c = cobj.get_alldata().squeeze() # get the concentration data from the file
        times = np.array(cobj.get_times()) # get the times and convert to array
    
        t_end_index = len(times)
        t_begin_index = int(0.9*len(times))
        climit = 1 # limit concentration, g/L
        
        time_break_lst = []
        rec_eff_lst = []
        cycle_n = np.arange(0, 10,1)
        c_arr = np.zeros((len(cycle_n)+1, ncol))
        c_store_all = np.zeros((len(cycle_n), nstepin+nstepout,ncol))
        c_arr[0] = np.array([cs] * ncol)
        c_prev = cs
        sim.simulation_data.verbosity_level = VerbosityLevel.quiet
        for index_cycle in cycle_n:
            # initial condition from prev time period
            gwt_ic = fp.mf6.ModflowGwtic(model=gwt, 
                                         strt=c_arr[index_cycle], # initial concentration
                                         ) 
            # here also change the injection after the first two years. 
        
            
            # write model, solve model, and read concentration data
            sim.write_simulation(silent=True)
            success, _ = sim.run_simulation(silent=True) 
            if success == 1:
                print(f'Model solved successfully for {index_cycle}', end="\r")
            else:
                print('Solve failed')
                break
            
            cobj = gwt.output.concentration() # get handle to binary concentration file
            c_i = cobj.get_alldata().squeeze() # get the concentration data from the file
            for itime in range(t_begin_index, t_end_index):
                if c_i[itime, 0] > climit:
                    time_break_lst.append(itime)
                    break
        
            c_arr[index_cycle+1] = c_i[itime - 1]
            c_store_all[index_cycle] = c_i
            rec_eff = ((times[itime - 1] - tin) * Q) / (tin * Q_in) # Q  needed as injection and extraction rates are not the same
            rec_eff_lst.append(rec_eff*100)

        
        fname = rf"figures/recovery_Q_cycles_kh{k}_npor{npor}.png"
        cycle_n_arr = np.array(cycle_n) + 1
        rec_eff_arr = np.array(rec_eff_lst)
        plt.figure()
        plt.plot(cycle_n_arr,rec_eff_arr*Q_tot/100,marker="o".label=f"\nkh={k}m/d, porosity={npor}")
        plt.axhline(Q_d,ls="--",color="C1",label="Design production volume")
        plt.ylabel("Recovery Volume [m^3]")
        plt.xlabel(r"Years after starting")
        plt.title(f"Recovery efficiency as cycle changes")
        plt.legend()
        plt.xticks(ticks=cycle_n_arr)
        plt.grid()
        plt.savefig(fname,bbox_inches="tight")
        n_run+=1
    plt.close()

print('DONE!')
