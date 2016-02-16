def chain(input_file, output_dir, potential, T, steps=None, rmsd=None, dt=0.004, rc=2.5):
    import rumd
    from rumdSimulation import rumdSimulation

    inp = input_file
    for Ti, di in zip(T, output_dir):
        for run in ['eq', 'av']:
            out_dir = "%s/%s" % (di, run)          
            if run == 'av':
                sim = single(inp, potential, T, dt, interval_energy=target_steps/1000, interval_config=target_steps/100)
                si = Simulation(sim, out_dir, target_steps=last_steps, restart=True)
            else:
                sim = single(inp, potential, T, dt, interval_energy=100000, interval_config=0)
                si = Simulation(sim, out_dir, target_rmsd=rmsd, target_steps=steps, restart=True)
            #from atooms.simulation import TargetRMSD
            #si.add(TargetRMSD(target_rmsd), Scheduler(period=100000))
            yield si

            # Chain with previous output file
            inp = '%s/final.xyz.gz' % out_dir
            last_steps = si.steps

def rx(input_file, output_dir_root, output_dir, potential, Tl, dt, nsteps, swap_period, config_period, seed, rmsd=None):
    """
    Tl: list of temperatures
    output_dir_rott: root directory for pt
    output_dir: list of output directories, should match lengths of Tl
    potential: function returning a list of rumd potentials
    """
    sim = rx(input_file, potential, Tl, dt)
    sim_adapter = [Simulation(s, output_dir_root + '/replica/%d' % i) for i, s in enumerate(sim)]
    rx_sim = ParallelTempering(output_dir_root, output_dir, zip(Tl, dt), sim_adapter, swap_period, seed=seed,
                               target_steps=nsteps, target_rmsd=rmsd, thermo_period=1, config_period=config_period, 
                               checkpoint_period=config_period, restart=True)
    yield rx_sim
