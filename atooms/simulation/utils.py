# TODO: this should be an API
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
