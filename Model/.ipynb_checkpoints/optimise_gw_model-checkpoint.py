# import the necessary packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
plt.rcParams['figure.figsize'] = (5, 3) # set default figure size
import flopy as fp  # import flopy and call it fp
import glob
import imageio
import IPython
import scipy
from tqdm import tqdm
import warnings

def run_gw_model(loss_factor, params,return_rec_eff=False):
    gwt, gwf, sim, Q_d, Q, d_injecting, cf, t_begin_index, t_end_index, ncol, nstepin,nstepout, cs, climit, tin =\
    params[0], params[1], params[2],params[3], params[4], params[5],params[6],params[7],params[8], params[9],params[10],params[11],params[12], params[13], params[14]
    """
    loss factors is the amount each year loses and thus what should be injected in that case.
        should be of length cycle_n and thus in this case 10 years.
    """
    try:
        assert len(loss_factor) == 10
    except AssertionError:
        print(loss_factors)
        raise AssertionError
    ### run model: 
    time_break_lst = []
    rec_eff_lst = []
    cycle_n = np.arange(0, 10,1)
    c_arr = np.zeros((len(cycle_n)+1, ncol))
    c_store_all = np.zeros((len(cycle_n), nstepin+nstepout,ncol))
    c_arr[0] = np.array([cs] * ncol)
    c_prev = cs
    
    lst_Q_in = []
    warnings.simplefilter("ignore") # warning in tqdm
    for index_cycle in tqdm(cycle_n):
        # initial condition from prev time period
        gwt_ic = fp.mf6.ModflowGwtic(model=gwt, 
                                     strt=c_arr[index_cycle], # initial concentration
                                     ) 
        Q_in = (Q_d * loss_factor[index_cycle]) / (d_injecting)# injection rate m^3/d
        lst_Q_in.append(Q_in)
        # here also change the injection after the first two years. 
        wellin = [[(0, 0, 0),  Q_in, cf]]   # [(layer, row, col), U, concentration]
        wellout = [[(0, 0, 0),  -Q, cf]] # specified concentration is not used, but must be specified 
        wel_spd = {0: wellin, 1: wellout} # stress period data for periods 0 and 1
        gwf_wel = fp.mf6.ModflowGwfwel(model=gwf, 
                               stress_period_data=wel_spd, 
                               auxiliary=['CONCENTRATION'],
                               pname='WEL1', # package name
                              )
        
        # write model, solve model, and read concentration data
        sim.write_simulation(silent=True)
        success, _ = sim.run_simulation(silent=True) 
        if success == 1:
            # print(f'Model solved successfully for {index_cycle}', end="\r")
            pass # replaced by tqdm
        else:
            print('Solve failed')
            break
        
        cobj = gwt.output.concentration() # get handle to binary concentration file
        c_i = cobj.get_alldata().squeeze() # get the concentration data from the file
        times = np.array(cobj.get_times()) # get the times and convert to array
        for itime in range(t_begin_index, t_end_index):
            if c_i[itime, 0] > climit:
                time_break_lst.append(itime)
                break
    
        c_arr[index_cycle+1] = c_i[itime - 1]
        c_store_all[index_cycle] = c_i
        rec_eff = ((times[itime - 1] - tin) * Q) / (tin * Q_in) # Q  needed as injection and extraction rates are not the same
        rec_eff_lst.append(rec_eff*100)

    Q_injected_total = np.array(lst_Q_in) * d_injecting
    rec_eff_arr = np.array(rec_eff_lst)
    recovered_volume = rec_eff_arr * Q_injected_total / 100
    volume_excess = recovered_volume - Q_d
    if return_rec_eff:
        return rec_eff_arr, volume_excess
    else:
        return volume_excess


def generate_plot(rec_eff_arr, volume_excess, params_plot):
    Q_d, cycle_n, loss_factor  = params_plot[0], params_plot[1], params_plot[2]
    cycle_n_arr = np.array(cycle_n) + 1
    fig, ax0 = plt.subplots(1,1)
    Q_injected_total = Q_d * loss_factor
    recovered_volume = rec_eff_arr * Q_injected_total / 100
    ax0.plot(cycle_n_arr,recovered_volume,marker=".",label="Model extraction volume")
    ax0.axhline(Q_d,ls="--",color="C1",label="Design production volume")
    ax0.set_ylabel("Volume [m^3]")
    ax0.set_xlabel(r"Time (years)")
    ax0.set_title(f"Injection vs extraction volumes over time")
    ax0.set_xticks(ticks=cycle_n_arr)
    ax0.grid()
    ax0.bar(cycle_n_arr,Q_injected_total,color="k",alpha=0.3,zorder=-1, label="Injection volume")
    ax0.legend()
    ax0.set_ylim(min(recovered_volume)-500,max(Q_injected_total)+500)
    # fig.savefig(fname,bbox_inches="tight")
    return fig