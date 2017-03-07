import rumd
potential = [rumd.Pot_LJ_12_6(cutoff_method = rumd.ShiftedPotential)]
potential[0].SetVerbose(False)
potential[0].SetParams(i=0, j=0, Epsilon=1.0, Sigma=1.0, Rcut=2.5)

