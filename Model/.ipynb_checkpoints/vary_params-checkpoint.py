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
from tqdm import tqdm
import warnings
warnings.simplefilter(action='ignore', category=DeprecationWarning)
import datetime
import xarray as xr

#### vary k & npor
# k_lst = [10, 20, 30, 40]
k_lst = [10, 35, 40]
# k_lst = [40]
# k_lst = [10]
# npor_lst = [0.2, 0.3, 0.4, 0.5]
npor_lst = [0.2,0.35, 0.5]
# npor_lst = [0.5]
n_runs = len(k_lst) * len(npor_lst)
alphaL = 0.5

#### vary alphaL
# npor = 0.35
# k = 30
# alphaL_lst = [0.5, 1, 1.5, 2]
# n_runs = len(alphaL_lst)

nlay = 20   #################
n_years = 3 ################
cycle_n = np.arange(0, n_years,1)
store_eff = np.zeros((n_runs,n_years))
n_run = 0
for k in tqdm(k_lst):
    for npor in npor_lst:
# for alphaL in tqdm(alphaL_lst):
        print(n_run,end="\n\n")
        Q_d = 40_000
        Q_tot = Q_d * 1.25 ## Change this later as we actually want to produce this, currentlty we still have 80% loss
        d_extrating = 62
        d_injecting = 365 - d_extrating ## change later for days with water excess
        # print(f'Need to pump {Q_tot/d_extrating:.2f}m^3/d to full fill demand') 
        
        # domain size and boundary conditions
        R = 200 # length of domain, m #############################################3333
        hR = 0 # head at r=R
        hL = hR

        # aquifer parameters
        # k = 40 # hydraulic conductivity, m/d -> ############ 10! 
        H = 20 # aquifer thickness, m - fixed
        # npor = 0.35 # porosity, generally 0.25 - 0.5  ############# 0.25
        
        # flow
        Q = Q_tot/d_extrating # extraction rate, m^3/d 
        Q_in = Q_tot / (d_injecting)# injection rate m^3/d
        
        # transport
        # alphaL = 0.5 # longitudinal dispersivity in horizontal direction, m - # something to check!! -> slides guest lecture
        alphaT = alphaL / 10 # transverse dispersivity in horizontal direction, m
        diffusion_coef = 0 # diffusion coeffcient
        
        # concentration
        cs = 30 # initial concentration, kg/m^3 (=g/L)
        cf = 0 # concentration injected water, kg/m^3 (=g/L)

        # buoyancy
        rhoref = 1000 # reference density, kg/m^3
        cref = 0 # reference concentration, kg/m^3
        drhodc = 0.7143  # Slope of the density-concentration line
        
        # space discretization
        delr = 0.2 # length of cell along row (in x-direction), m
        delc = 1 # width of cells normal to plane of flow (in y-direction), m
        z = np.linspace(0, -H, nlay + 1) # top and bottom(s) of layers
        zc = 0.5 * (z[:-1] + z[1:]) # center so cells, used for contouring
        nrow = 1 # number of rows
        ncol = round(R / delr) # number of columns
        
        # radialize parameters:
        theta = 2 * np.pi
        r = np.cumsum(delr * np.ones(ncol)) - 0.5 * delr # rightside of cell minus half the cell length
        krad = k * r * theta * np.ones((nlay, nrow, ncol))
        nporrad = npor * r * theta * np.ones((nlay, nrow, ncol))
        
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
        modelname = f'modelrad__{k}_{str(npor).split(".")[-1]}' # name of model
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
        ################# GWF #######################
        # groundwater flow model
        gwf = fp.mf6.ModflowGwf(simulation=sim, # add to simulation called sim
                                modelname=gwfname, # name of gwf model
                                save_flows=True, # make sure all flows are stored in binary output file
                               )
        
        # iterative model solver
        gwf_ims  = fp.mf6.ModflowIms(simulation=sim, # add to simulation called sim
                                     filename=gwf.name + '.ims', # file name to store ims
                                     linear_acceleration="BICGSTAB", # use BIConjuGantGradientSTABalized method
                                     inner_dvclose=1e-6
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
                                        k33=krad / 10,
                                        save_flows=True, # save the flow for all cells
                                       )
            
        # initial condition
        gwf_ic = fp.mf6.ModflowGwfic(model=gwf, 
                                     strt=hR, # initial head used for iterative solution
                                    )
        
        # wells
        wellin = []
        wellout = []
        for ilay in range(nlay):
            wellin.append([(ilay, 0, 0),  Q / nlay, cf])  # [(layer, row, col), U, concentration]
            wellout.append([(ilay, 0, 0),  -Q / nlay, cf]) # specified concentration is not used, but must be specified 
        wel_spd = {0: wellin, 1: wellout} # stress period data for periods 0 and 1
        gwf_wel = fp.mf6.ModflowGwfwel(model=gwf, 
                                       stress_period_data=wel_spd, 
                                       auxiliary=['CONCENTRATION'],
                                       pname='WEL1', # package name
                                      )
        
        # constant head 
        chd = []
        for ilay in range(nlay):
            chd.append([(ilay,  0,  ncol-1), hL, cs]) # [(layer, row, col), head, concentration]
        chd_spd  = {0: chd, 1: chd}    # Stress period data
        gwf_chd = fp.mf6.ModflowGwfchd(model=gwf, 
                                       stress_period_data=chd_spd, 
                                       auxiliary=['CONCENTRATION'],
                                       pname='CHD1', # package name
                                      )
        
        # buoyancy
        buy = fp.mf6.ModflowGwfbuy(model=gwf,
                                   packagedata=[0, drhodc, cref, gwtname, 'CONCENTRATION'], # [conc 1 species - 0= salt, drhodc, cref, gwtname, name]
                                   denseref=rhoref, # reference concentration
                                   nrhospecies=1, # number of species
                                   density_filerecord=f"{gwf.name}.dst", # file name
                                   pname='BUY1', 
                                  )
            
        # output control
        oc = fp.mf6.ModflowGwfoc(model=gwf, 
                                 saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")], # what to save
                                 budget_filerecord=f"{gwfname}.cbc", # file name where all budget output is stored
                                 head_filerecord=f"{gwfname}.hds", # file name where all head output is stored
                                )

        ################# GWT #######################
        # groundwater transport model
        gwt = fp.mf6.ModflowGwt(simulation=sim, 
                                modelname=gwtname, # name of groundwater transport model
                               )
        
        # iterative model solver
        gwt_ims  = fp.mf6.ModflowIms(simulation=sim,
                                     filename=gwt.name + '.ims', # must be different than file name of gwf model ims
                                     linear_acceleration="BICGSTAB",
                                     inner_dvclose=1e-6
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
        
        # sim.write_simulation(silent=True)
        # success, _ = sim.run_simulation(silent=True) 
        # if success == 1:
        #     print('Model solved successfully',end=" ")
        # else:
        #     print('Solve failed')
        
        # cobj = gwt.output.concentration() # get handle to binary concentration file
        # c = cobj.get_alldata().squeeze() # get the concentration data from the file
        # times = np.array(cobj.get_times()) # get the times and convert to array
    
        t_end_index = nstepin + nstepout
        t_begin_index = nstepin
        climit = 1 # limit concentration, g/L

        time_break_lst = []
        rec_eff_lst = []
        cycle_n = np.arange(0, n_years,1)
        c_arr = np.zeros((len(cycle_n)+1, nlay, ncol))
        c_store_all = np.zeros((len(cycle_n), nstepin+nstepout,nlay,ncol))
        c_arr[0] = np.ones((nlay,ncol)) * cs
        c_prev = cs
    
        now = datetime.datetime.now()
        print(f'Start {nlay} layers at {now.hour}:{now.minute}')
        end = now + datetime.timedelta(minutes=(nlay-5)*n_years)
        print(f'Expected runtime end at {end.hour}:{end.minute}')
        for index_cycle in tqdm(cycle_n):
            # initial condition from prev time period
            gwt_ic = fp.mf6.ModflowGwtic(model=gwt, 
                                         strt=c_arr[index_cycle], # initial concentration
                                         ) 
            # here also change the injection after the first two years. 
        
            
            # write model, solve model, and read concentration data
            sim.write_simulation(silent=True)
            success, _ = sim.run_simulation(silent=True) 
            if success == 1:
                # print(f'Model solved successfully for {index_cycle}', end="\r")
                pass
            else:
                print('Solve failed')
                break
            
            cobj = gwt.output.concentration() # get handle to binary concentration file
            c_i = cobj.get_alldata().squeeze() # get the concentration data from the file
            for itime in range(t_begin_index, t_end_index):
                if c_i[itime,:, 0].mean() > climit:
                    time_break_lst.append(itime)
                    break
        
            c_arr[index_cycle+1] = c_i[itime - 1,:]
            c_store_all[index_cycle] = c_i
            rec_eff = ((times[itime - 1] - tin) * Q) / (tin * Q_in) # Q  needed as injection and extraction rates are not the same
            rec_eff_lst.append(rec_eff*100)
            print(f'{k}_{npor}_{index_cycle}_{itime}_{rec_eff}')

        print(rec_eff_lst)
        print(time_break_lst)
        cycle_n_arr = np.array(cycle_n) + 1
        rec_eff_arr = np.array(rec_eff_lst)
        print(rec_eff_arr)
        store_eff[n_run,:] = rec_eff_arr
        n_run+=1
        time = str(datetime.datetime.now())[:-10].replace(":","_")
        fname = fr'output/store_concentrations_k-{k}_npor-{npor}_alphaL-{alphaL}-nlay-{nlay}_nyears-{nyears}_{time}.nc'
        ds = xr.DataArray(c_store_all,dims=['year','tstep','layer','r'])
        ds.to_netcdf(fname,engine="netcdf4")

time = str(datetime.datetime.now())[:-10].replace(":","_")
np.savetxt(fr'output/model_output_{time}.txt',store_eff,delimiter=",")
params = np.array(k_lst + npor_lst)
np.savetxt(fr'output/parameter_output_{time}.txt',params,delimiter=",")
print('DONE!')
