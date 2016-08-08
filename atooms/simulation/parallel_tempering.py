# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

import os
import sys
import math
import random
import numpy
import datetime

from atooms.simulation import Simulation, WriterCheckpoint
from atooms.simulation import log
from atooms.utils import rank, size, comm, barrier
from atooms.utils import rmd, rmf, mkdir

class WriterConfig(object):

    def __str__(self):
        return 'config'

    def __call__(self, e):
        log.debug('write config %d' % e.steps)
        for i in e.my_replica:
            irx = e.state[i]
            # If mute_config_except is None or it contains the state we write it
            # otherwise we do not write anything. 
            # Ex: only write the lowest T state -> mute_config_except=[0]
            #     write all states -> mute_config_except=None (default)            
            if e.mute_config_except is None or \
               irx in e.mute_config_except:
                # TODO: this is what we should do if we wanted to use the underlying backend of the replica with state irx
                # e.sim[irx].writer_config(e.sim[irx])
                # Issue is how to avoid writing velocities
                # Issue is we should make sure the backend itself doesn't call the writer...
                # One way to do that is to delegate to backend entirely, i.e. pt should not add writer_config
                # as observer. We'd have to make sure that simulation backends write configurations at the right
                # interval... this will be against the idea that the backend is just there to propagate the dynamics!
                #
                # So we need two things here: a trajectory class that is appropriate to the backend
                # (actually, any trajectory would be fine as long as the backend implements System interface)
                # and an output file (or output directory). In writer_thermo here this is set by the writer
                # Perhaps we should define it here. We could grab the output_path from pt instance and the suffix from
                # the trajectory class
                #
                # The trajectory class is taken from the simulation backend
                trj_cls = e.sim[irx].trajectory
                d = e.dir_state_out[irx]
                with trj_cls(os.path.join(d, 'trajectory.' + trj_cls.suffix), 'a') as t:
                    t.exclude(['velocity'])
                    t.write(e.replica[i], e.steps)

class WriterCheckpointPT(object):

    def __str__(self):
        return 'checkpoint'

    def __call__(self, e):
        e.write_checkpoint()
        
class WriterThermo(object):

    def __str__(self):
        return 'thermo'

    def __call__(self, e):
        # Since we grab steps from simulations, we must gather them first
        # We could have each process write down its replicas and make it more efficient, see write_state()
        u = numpy.ndarray(e.nr)
        u_l = numpy.ndarray(len(e.my_replica))
        k = numpy.ndarray(e.nr)
        k_l = numpy.ndarray(len(e.my_replica))
        rmsd = numpy.ndarray(e.nr)
        rmsd_l = numpy.ndarray(len(e.my_replica))
        steps = numpy.ndarray(e.nr, dtype=int)
        steps_l = numpy.ndarray(len(e.my_replica), dtype=int)
        for i, ri in enumerate(e.my_replica):
            u_l[i] = e.replica[ri].potential_energy()
            k_l[i] = e.replica[ri].kinetic_energy()
            rmsd_l[i] = e.sim[ri].rmsd
            steps_l[i] = e.sim[ri].steps

        if size > 1:
            comm.Gather(u_l, u, 0)
            comm.Gather(k_l, k, 0)
            comm.Gather(rmsd_l, rmsd, 0)
            comm.Gather(steps_l, steps, 0)
        else:
            u = u_l
            k = k_l
            rmsd = rmsd_l
            steps = steps_l

        if rank == 0:
            #e.write_replica(rmsd, steps)
            #e.write_state(u, k, steps)

            # Loop over replicas
            for i in range(e.nr):
                f = e.output_path + '/replica/%d.out' % i
                with open(f, 'a') as fh:
                    # In which state is physical replica i ?
                    fh.write('%d %d %d %g\n' % (e.steps, steps[i], e.state[i], rmsd[i]))

            # Loop over states
            for i in range(e.nr):
                f = e.output_path + '/state/%d.out' % i
                with open(f, 'a') as fh:
                    # Which replica is in state i? What is its energy?
                    irep = e.replica_id[i]
                    fh.write('%d %d %d %.6g %.6g\n' % (e.steps, steps[i], irep, 
                                                       u[irep], k[irep]))


class StateTemperature(object):

    """Callback to set the state of a simulation, once PT has swapped states"""

    def __call__(self, system, params):
        # Note that when swapping certan thermostats, the internal state should be reset.
        # It is OK for RUMD to only set the temperatures from the params list
        try:
            system.thermostat.temperature = params[0]
            system.dynamics.timestep = params[1]
        except TypeError:
            system.thermostat.temperature = params


# Design:

# Trajectory are defined internally by Simulation
# The actual backend might be changed by the user at run time

# Trajectory are usually one per file but rumd uses a directory based logic
# From a file name we can always reconstruct the dirname,
# the other way around we need a hint how to construct the path, like a base
# From this discussion it looks like using file name in constructor is always best ... (even for sim)

# It's not the replica simulation responsibility to handle the trajectory backend
# Each simulation has its own trajectory backend, which can be changed of course

# In RX we consider a step=block since the number of steps may depend on the replica

class ParallelTempering(Simulation):

    """
    Parallel tempering simulation.
    """

    _WRITER_THERMO = WriterThermo
    _WRITER_CONFIG = WriterConfig    
    _WRITER_CHECKPOINT = WriterCheckpointPT

    # TODO: can we pass a NamedTuple like params (perhaps coming straight from argparse) instead? It would be more compact, more extensible, but less safe (we should check for missing variables and sane defaults)
    def __init__(self, sim, params, output_path,
                 swap_interval=0, seed=10, update=StateTemperature, fmt='T%.4f',
                 mute_config_except=None,
                 steps=None, rmsd=None,
                 thermo_interval=0, thermo_number=0, 
                 config_interval=0, config_number=0,
                 checkpoint_interval=0, checkpoint_number=0,
                 restart=False, dryrun=False, swap_scheme='alternate_all'):
        # Note we ignore the backend here
        Simulation.__init__(self, output_path=output_path,
                            steps=steps, rmsd=rmsd,
                            thermo_interval=thermo_interval, thermo_number=thermo_number, 
                            config_interval=config_interval, config_number=config_number,
                            checkpoint_interval=checkpoint_interval, checkpoint_number=checkpoint_number,
                            restart=restart)
        self.params = params
        self.sim = sim
        self.targeter_rmsd_period = 10
        # TODO: drop variables, make acceptance a callback
        self.variables = ['T']
        self.update = update()
        self.steps_block = swap_interval
        self.mute_config_except = mute_config_except
        self.nr = len(params)
        self.seed = seed
        self.dryrun = dryrun
        self.swap_scheme = swap_scheme
        random.seed(self.seed)

        # Get physical replicas (systems) from simulation instances.
        # These are references: they'll follow the simulations 
        self.replica = [s.system for s in sim]
        self.system = self.replica[0]

        # Sanity check
        if not (self.nr == len(sim)):
            raise ParameterError('n. of backends must match n. of states (%d, %d)' % (self.nr, len(sim)))

        # Define output files
        self.file_log = self.output_path + '/pt.log'
        # For each thermodynamic state, info on the replica which has it
        self.file_state_out = [self.output_path + '/state/%d.out' % i for i in range(self.nr)]
        # For each physical replica, info on the state in which it is
        self.file_replica_out = [self.output_path + '/replica/%d.out' % i for i in range(self.nr)]
        # This will be set as output_path of the backend in a moment
        if fmt is None:
            self.dir_state_out = [self.output_path + '/state/%d' % i for i in range(self.nr)]
        else:
            self.dir_state_out = [self.output_path + ('/state/' + fmt) % p for p in self.params]
        for s, d in zip(self.sim, self.dir_state_out):
            # This variable should be None (or we are overwriting something)
            s.output_path = d

        self.replica_id = range(self.nr)
        # TODO: remember to checkpoint them
        self.accepted = [0.0 for i in range(self.nr)]
        self.attempts = [0.0 for i in range(self.nr)]
        # This is used to swap odd/even replicas in turn
        self.offset = 0

        # Distribute physical replicas in parallel.
        # Each process gets a bunch of replicas to evolve
        # replica_id contains their state id's.
        # We could as well distributes states, which would allow
        # for other optimizations.
        np = self.nr / size
        ni = rank * np
        nf = (rank+1) * np
        self.my_replica= range(ni, nf)
        self.process_with_replica = numpy.array(range(self.nr))
        for irank in range(size):
            ni = irank * np
            nf = (irank+1) * np
            for nr in range(ni, nf):
                self.process_with_replica[nr] = irank
        barrier()

        # Listify steps per block.
        # Note that this is per state. If we set it in the sim instances
        # We should update them all the time.
        if not type(swap_interval) is list:
            self.steps_block = [swap_interval] * self.nr
        # If we provided 1 block interval in list, upgrade it to nr length
        if len(self.steps_block) != self.nr:
            if len(self.steps_block) == 1:
                self.steps_block = swap_interval * self.nr
            else:
                raise ValueError('number of blocks does not match n. replica')

    @property
    def rmsd(self):
        """In parallel tempering simulation we define rmsd as the minimum one"""
        # RMSD must be known on all processes, thus an allgather is needed
        # This might be optimized by targeting rmsd dynamically.
        rmsd_l = numpy.ndarray(len(self.my_replica))
        for i, ri in enumerate(self.my_replica):
            rmsd_l[i] = self.sim[ri].rmsd
        if size > 1:
            rmsd = numpy.ndarray(self.nr)
            comm.Allgather(rmsd_l, rmsd)
            return min(rmsd)
        else:
            return min(rmsd_l)

    def clean_files(self):
        for f in \
                self.file_state_out + \
                self.file_replica_out:
            rmf(f)

    def write_log(self):
        """ Dump state parameters """
        with open(self.file_log, 'w') as fh:
            for i, T in enumerate(self.params):
                try:
                    fmt = ' %s' * len(T)
                    fh.write(('%d' + fmt + '\n') % tuple([i] + list(T)))
                except:
                    fh.write('%d %s\n' % (i, T))

    # TODO: drop these two as they are part of WriterThermo
    def write_replica(self, msd, step):
        """ Dump output info on a physical replica """
        # Loop over replicas
        for i in range(self.nr):
            with open(self.file_replica_out[i], 'a') as fh:
                # In which state is physical replica i ?
                fh.write('%d %d %d %g\n' % (self.steps, step[i], self.state[i], msd[i]))

    def write_state(self, u, k, step):
        """ Dump output info on a thermodynamic state """
        # TODO: we could write state atomically, which would allow parallelization and remove communications
        # Loop over states       
        log.debug('rx step=%s replicas(state)=%s' % (step[0], self.replica_id))

        for i in range(self.nr):
            with open(self.file_state_out[i], 'a') as fh:
                # Which replica is in state i? What is its energy?
                fh.write('%d %d %d %g %g\n' % (self.steps, step[i], self.replica_id[i], 
                                               u[self.replica_id[i]], k[self.replica_id[i]] )) #, self.acceptance(i)))

    # def read_checkpoint(self):
    #     """ Checkpoint """
    #     for i in self.my_replica:
    #         f = self.file_state_out[i] + '.chk'
    #         with open(f, 'r') as fh:
    #             self.state[i] = fh.read()
    #         self.replica[i], steps = self.trajectory[i].read_checkpoint()

    def write_checkpoint(self):
        """Checkpoint replicas via simulation backends as well as the
        thermodynamic states in which the replicas found themselves.
        """
        # TODO: make this more robust. Only if all check points (state and replica) have been written we can safely restart!
        # We should therefore first keep the old checkpoints, create new files and then move (which is quick).
        # Additionally we should check consistency upon reading (i.e. all checkpoints should belong to the same step) and fail otherwise
        log.debug('write checkpoint %d' % self.steps)
        # Note: offset and step are redundant, since they are global
        for i in self.my_replica:
            with open(self.file_replica_out[i] + '.chk', 'w') as fh:
                fh.write('%d\n' % self.state[i])
                fh.write('%d\n' % self.steps)
                fh.write('%d\n' % self.offset)
            self.sim[i].write_checkpoint()

    def check(self):
        for i in range(self.nr):
            T = float(self.replica[i].thermostat.temperature)
            if abs(self.params[self.state[i]] - T) > 1e-5:
                log.error('replica %d state %d at T=%s has thermostat at %s, delta %f' % (i, self.state[i], self.params[self.state[i]], T, abs(self.params[self.state[i]] - T)))
                raise RuntimeError

    def run_pre(self):
        Simulation.run_pre(self)

        # If we do not restart, we clear up everything in the base
        if not self.restart:
            rmf(self.output_path + '/pt.log')
            rmd(self.output_path + '/state')
            rmd(self.output_path + '/replica')

        # Make sure base directories exist
        mkdir(self.output_path)
        mkdir(self.output_path + '/state')
        mkdir(self.output_path + '/replica')
        # Make sure output directories exist, even if we dont write config (for checkpoint)
        for d in self.dir_state_out:
            mkdir(d)

        if self.restart:
            # TODO: steps should all be equal, we should check
            # This must be done by everybody. Otherwise, each process
            # should read its replicas and then gather.
            for i in range(self.nr):
                self.sim[i].restart = True
                # This is basically a read_checkpoint()
                f = self.file_replica_out[i] + '.chk'
                if os.path.exists(f):
                    with open(f, 'r') as fh:
                        istate = int(fh.readline())
                        self.replica_id[istate] = i
                        self.steps = int(fh.readline())
                        self.offset = int(fh.readline())
                    log.debug('pt restarting replica %d at state %d from step %d' % 
                                  (i, istate, self.steps))
        # Restarting is handled by the simulation instance.
        # The only glitch for now is that the checkpoint file ends up
        # in the directory corresponding to the initial state of the replica
        for i in self.my_replica:
            self.sim[i].run_pre()

        self.write_log()

    def __str__(self):
        return 'Parallel tempering simulation'

    def wall_time_per_step(self):
        """Return the wall time in seconds per step (over all replicas)"""
        return self.elapsed_wall_time() / (self.steps-self.initial_steps) / sum(self.steps_block)

    def _report(self):
        log.info('backend: %s' % self.sim[0])
        log.info('output path: %s' % self.output_path)
        log.info('number of replicas: %d' % self.nr)
        log.info('number of processes: %d' % size)
        if self.dryrun:
            log.info('** this is a dry run (no actual swaps) **')
        if self.steps_block[0] == self.steps_block[-1]:
            log.info('swap interval: %d' % self.steps_block[0])
        barrier()
        log.info('process %s has replicas: %s at state %s' % (rank, self.my_replica, [self.state[i] for i in self.my_replica]), extra={'rank':'all'})

    def _report_end(self):
#        log.info('final minimum acceptance: %.2f' % min([self.acceptance(i) for i in range(self.nr)]))
        log.info('final minimum rmsd: %.2f' % self.rmsd)
        log.info('wall time [s]: %.1f' % self.elapsed_wall_time())
        log.info('average TSP [s/step/particle]: %.2e' % (self.wall_time_per_step_particle()))
        log.info('simulation ended on: %s' % datetime.datetime.now().strftime('%Y-%m-%d at %H:%M'))

    def run_until(self, nsteps):
        """Run until nsteps, which is the number of PT steps, i.e. a block of several steps"""
        # In general, nsteps will be 1, but with dryrun we may do more steps at once
        # Evolve over my physical replicas.
        log.debug('run until %d' % nsteps)
        for i in self.my_replica:
            # This will evolve physical replica[i] for the number
            # of steps prescribed for its state times the number of blocks 
            # (relevant only for dryrun)
            n = nsteps * self.steps_block[self.state[i]]
            #n = self.sim[i].steps + (nsteps-self.steps) * self.steps_block[self.state[i]]
            # TODO: rather use an additional level (debug_all, info_all) than extra dict
            log.debug('evolve replica %d on GPU %d for %d until %d' % (i, rank, nsteps,
                                                                       n), extra={'rank':'all'})
            log.debug('replica %d state %d formally at T=%s has thermostat T %s' % \
                      (i, self.state[i], self.params[self.state[i]], self.replica[i].thermostat.temperature))
            self.sim[i].run_until(n)
        # Its up to us to update our steps
        self.steps = nsteps
        # Attempt to exchange replicas.
        if not self.dryrun:
            self.exchange(self.replica)

    def run_end(self):
        # TODO: not called anymore
        for i in self.my_replica:
            self.sim[i].run_end()
        barrier()

    @property
    def state(self):
        # TODO: make state a private list. Done this way is unsafe
        # because python will not raise an attribute error if we try
        # to modify an item of the list!
        s = [0] * self.nr
        for i in range(self.nr):
            s[self.replica_id[i]] = i
        return s

    def acceptance(self, state):
        if self.attempts[state] > 0:
            return self.accepted[state] / self.attempts[state]
        else:
            return 0.0

    def __repr__(self, state):
        return 'state %d <-> %d, replicas %d <-> %d, [%d/%d], prob %.4f]' % \
            (state, state+1, self.replica_id[state], self.replica_id[state+1], self.accepted[state], self.attempts[state], self.__prob)

    def _swap(self, system, state, nn_state, r_i=None, r_j=None):
        # TODO: ri and rj unused
        tmp = self.replica_id[state]
        self.replica_id[state] = self.replica_id[nn_state]
        self.replica_id[nn_state] = tmp

    def _process_with(self, state):
        return self.process_with_replica[self.replica_id[state]]

    def exchange(self, system):
        log.debug("exchange rx")
        if self.variables == ['T']:
            if self.swap_scheme == 'random_one':
                self.__exchange_T_one(system)
            elif self.swap_scheme == 'alternate_all':
                self.__exchange_T_parallel(system)
        else:
            raise ValueError(self.variables)

    def __exchange_T_parallel(self, system):
        # Loop over ALL states, find out which states are sending out messages
        barrier()
        sender = []
        focused = []
        for state in range(self.nr):
            if self.replica_id[state] in self.my_replica:
                focused.append(state)
            else:
                continue
            if (state + self.offset) % 2 == 0:
                sender.append(state)
        
        for state in range(self.nr):
            if not state in focused:
                continue
            if state in sender:
                if state < self.nr-1:
                    self.__exchange_T_comm(state, state+1, system)
            else:
                if state > 0:
                    self.__exchange_T_comm(state, state-1, system)

        # Update offset and replica ids
        self.offset = (self.offset+1) % 2

        # Once states have changed, we must update replica id's across processes
        state_tmp = numpy.array(range(self.nr))
        state_new = numpy.array([self.state[r] for r in self.my_replica])
        if size > 1:
            comm.Allgather(state_new, state_tmp)
        else:
            state_tmp = state_new

        # After setting replica_id we can safely use self.state
        for i, s in enumerate(state_tmp):
            self.replica_id[s] = i

        # -------------------------------
        # Update state of replicas
        for i, s in enumerate(self.state):
            self.update(system[i], self.params[s])
        # -------------------------------

    def __exchange_T_one(self, system):
        """Swap attempt one replica at a time"""
        # Get a random state i between 0 and nr-2
        barrier()
        if rank == 0:
            state = random.randint(0, self.nr-2)
        if size > 1:
            comm.Broadcast(state)

        # Attempt swap between states i and i+1, only concerned
        # processes call the exchange function
        if self.replica_id[state] in self.my_replica or \
           self.replica_id[state+1] in self.my_replica:
            if self._swap_attempt(system, state, state+1):
                self._swap(system, state, state+1)
                # TODO: fix accepted counter

        # TODO: no need to update anything if swap was unsuccesfull
        # Once states have changed, we must update replica id's across processes
        state_tmp = numpy.array(range(self.nr))
        state_new = numpy.array([self.state[r] for r in self.my_replica])
        if size > 1:
            comm.Allgather(state_new, state_tmp)
        else:
            state_tmp = state_new

        # After setting replica_id we can safely use self.state
        for i, s in enumerate(state_tmp):
            self.replica_id[s] = i

        # Update state of replicas using update() function
        for i, s in enumerate(self.state):
            self.update(system[i], self.params[s])

    def _swap_attempt(self, system, state_i, state_j):
        """This function should be called by processes that hold replicas at i and j"""
        # TODO: add attempts
        # Compute potential energy
        r_i = self.replica_id[state_i]
        r_j = self.replica_id[state_j] # This should be up to date
        u_i = numpy.array([system[r_i].potential_energy()])

        # Get processes holding replicas with states i and j
        rank_i = self._process_with(state_i)
        rank_j = self._process_with(state_j)

        if rank_i == rank_i == rank:
            log.debug('comm rank %d on same rank %d -> %d' % (rank, state_i, state_j))
            u_j = numpy.array([system[r_j].potential_energy()])
            ran = numpy.array([random.random()])
        else:
            # TODO: use sendrcv
            if rank == rank_i:
                log.debug('comm rank %d on /= rank %d -> %d' % (rank, state_i, state_j))
                # I am on the left, I send first.
                u_j = numpy.array([0.0])
                comm.Send(u_i, nn, 10)
                comm.Recv(u_j, nn, 11)
                ran = numpy.array([random.random()])
                comm.Send(ran, nn, 12)
            elif rank == recv_rank:
                log.debug('comm rank %d on /= rank %d <- %d' % (rank, state_j, state_i))
                # I am on the right I receive first
                u_j = numpy.array([0.0])
                comm.Recv(u_j, nn, 10)
                comm.Send(u_i, nn, 11)
                ran = numpy.array([0.0])
                comm.Recv(ran, nn, 12)
            else:
                raise ValueError("rank %d should not call this fct" % (rank))

        # Temperatures
        # TODO: encapsulate in exchange attempt callback
        try:
            T_i = self.params[state_i][0]
            T_j = self.params[state_j][0]
        except:
            T_i = self.params[state_i]
            T_j = self.params[state_j]

        # Store current probability term
        log.debug("comm rank %d uj=%g uu=%g Ti=%g Tj=%g" % (rank,  u_j[0], u_i[0], T_i, T_j))
        try:
            self.__prob = math.exp(-(u_j[0]-u_i[0])*(1/T_i-1/T_j))
        except:
            log.error('acceptance test failed uj=%g uj=%g Ti=%g Tj=%g' % (u_j[0], u_i[0], T_i, T_j))
            raise

        # Test if we can swap states of replicas
        log.debug("comm rank %d sawpping ran %g prob %g => %s" % (rank, ran[0], self.__prob, ran[0]<self.__prob))
        return ran[0] < self.__prob
            

    def __exchange_T_comm(self, my_state, nn_state, system):
        # TODO: add attempts
        exchange_attempt = False
        me = rank 
        # Compute potential energy
        r_i = self.replica_id[my_state]
        r_j = self.replica_id[nn_state] # This should be up to date
        u_i = numpy.array([system[r_i].potential_energy()])
        # Get process that has my nearest neighboring state
        nn = self._process_with(nn_state)

        if nn == rank:
            if my_state > nn_state:
                return
            log.debug('comm rank %d on same rank %d -> %d' % (rank, my_state, nn_state))
            u_j = numpy.array([system[r_j].potential_energy()])
            ran = numpy.array([random.random()])
        else:
            if my_state < nn_state:
                log.debug('comm rank %d on /= rank %d -> %d' % (rank, my_state, nn_state))
                # I am on the left, I send first.
                u_j = numpy.array([0.0])
                comm.Send(u_i, nn, 10)
                comm.Recv(u_j, nn, 11)
                ran = numpy.array([random.random()])
                comm.Send(ran, nn, 12)
            elif my_state > nn_state:
                log.debug('comm rank %d on /= rank %d <- %d' % (rank, my_state, nn_state))
                # I am on the right I receive first
                u_j = numpy.array([0.0])
                comm.Recv(u_j, nn, 10)
                comm.Send(u_i, nn, 11)
                ran = numpy.array([0.0])
                comm.Recv(ran, nn, 12)
            else:
                raise ValueError("my_state %d == nn_state %d" % (my_state, nn_state))

        # Temperatures
        # TODO: encapsulate in exchange attempt callback
        try:
            T_i = self.params[my_state][0]
            T_j = self.params[nn_state][0]
        except:
            T_i = self.params[my_state]
            T_j = self.params[nn_state]
        # Store current probability term
        log.debug("comm rank %d uj=%g uu=%g Ti=%g Tj=%g" % (rank,  u_j[0], u_i[0], T_i, T_j))
        try:
            self.__prob = math.exp(-(u_j[0]-u_i[0])*(1/T_i-1/T_j))
        except:
            log.error('acceptance test failed uj=%g uj=%g Ti=%g Tj=%g' % (u_j[0], u_i[0], T_i, T_j))
            raise
        # Test if we can swap states of replicas
        log.debug("comm rank %d sawpping ran %g prob %g => %s" % (rank,  ran[0], self.__prob, ran[0]<self.__prob,))
        if (ran[0] < self.__prob):
            self._swap(system, my_state, nn_state, r_i, r_j)
            # TODO: fix accepted counter

    def __exchange_T(self, system):
        """Simple non-parallel version"""
        # Sanity check
        if not len(system) == self.nr:
            raise ValueError

        # Loop over states defined by self.params with a jump of 2
        # The offset to swap only odd / even states at each attempt 
        for state in range(self.offset, self.nr-1, 2):
            self.attempts[state] += 1
            self.attempts[state+1] += 1
            T_i = self.params[state]
            T_j = self.params[state+1]
            # Index of physical replicas having states to swap
            r_i = self.replica_id[state]
            r_j = self.replica_id[state+1]
            # Get energies of replicas
            u_i = system[r_i].potential_energy()
            u_j = system[r_j].potential_energy()
            # Store current probability term
            self.__prob = math.exp(-(u_j-u_i)*(1/T_i-1/T_j))
            ran = random.random()
            # Test if we can swap states of replicas
            if (ran < self.__prob):
                self._swap(system, state, r_i, r_j)
                self.accepted[state] += 1
                self.accepted[state+1] += 1

        # Update offset
        self.offset = (self.offset+1) % 2



