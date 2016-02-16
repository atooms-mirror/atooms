def single(input_file, potential, T, dt, interval_energy="none", interval_config="none"):
    import rumd
    from rumdSimulation import rumdSimulation

    sim = rumdSimulation(input_file)
    for pot in potential():
        sim.AddPotential(pot)

    itg = rumd.IntegratorNVT(targetTemperature=T, timeStep=dt)
    sim.SetIntegrator(itg)
    sim.SetMomentumResetInterval(10000)

    if period is not None:
        sim.SetOutputScheduling("energies","linear",interval=interval_energy)
        sim.SetOutputScheduling("trajectory","linear",interval=interval_config)

    # Output dir??

    return sim
    #yield si

def rx(input_file, potential, Tl, dt):
    """
    Tl: list of temperatures
    output_dir_rott: root directory for pt
    output_dir: list of output directories, should match lengths of Tl
    potential: function returning a list of rumd potentials
    """
    import rumd
    from atooms.utils import size, rank, barrier
    from rumdSimulation import rumdSimulation
    
    # Create simulation and integrators
    for i in range(size):
        if i == rank:
            sim = [rumdSimulation(f) for f in input_file]
        barrier()
    igt = [rumd.IntegratorNVT(targetTemperature=T, timeStep=dti) for T, dti in zip(Tl, dt)]

    # Add potentials
    for s in sim:
        s.SetOutputScheduling("energies", "none")
        s.SetOutputScheduling("trajectory", "none")
        for p in potential():
            s.AddPotential(p)

    for s, i in zip(sim, igt):
        s.SetMomentumResetInterval(0)
        s.SetIntegrator(i)

    return sim
