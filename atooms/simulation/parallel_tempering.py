# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

import os
import sys
import math
import random
import numpy
import logging

from atooms.simulation import Simulation
from atooms.utils import rmd, rmf, mkdir

log = logging.getLogger()

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    logging.info('found mpi4py %d %d' % (rank, size))
except:
    rank = 0
    size = 1
    logging.info('mpi4py not found')

def barrier():
    if size > 1:
        comm.barrier()

class WriterConfig(object):

    def __call__(self, e):
        logging.debug('writer config')
        for i in e.my_replica:
            irx = e.state[i]
            # If the output directory for state irx is None
            # we do not write configurations
            if e.output_path[irx]:
                with e.trajectory(e.output_path[irx]+'/'+e.sim[irx].base_output, 'a') as t:
                    t.write_sample(e.replica[i], e.steps, ignore=['vel'])

class WriterCheckpoint(object):

    def __call__(self, e):
        e.write_checkpoint()
        
class WriterThermo(object):

    def __call__(self, e):
        logging.debug('writer thermo')

        # Since we grab steps from simulations, we must gather them first
        # We could have each process write down its replicas and make it more efficient, see write_state()
        u = numpy.ndarray(e.nr)
        u_l = numpy.ndarray(len(e.my_replica))
        rmsd = numpy.ndarray(e.nr)
        rmsd_l = numpy.ndarray(len(e.my_replica))
        steps = numpy.ndarray(e.nr, dtype=int)
        steps_l = numpy.ndarray(len(e.my_replica), dtype=int)
        for i, ri in enumerate(e.my_replica):
            u_l[i] = e.replica[ri].potential_energy()
            rmsd_l[i] = e.sim[ri].rmsd
            steps_l[i] = e.sim[ri].steps

        if size > 1:
            comm.Gather(u_l, u, 0)
            comm.Gather(rmsd_l, rmsd, 0)
            comm.Gather(steps_l, steps, 0)
        else:
            u = u_l
            rmsd = rmsd_l
            steps = steps_l

        if rank == 0:
            e.write_replica(rmsd, steps)
            e.write_state(u, steps)

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
    
    It relies on ReplicaExchange for handling replicas and swaps and on
    a list of Simulation objects *sim*.
    """

    _WRITER_THERMO = WriterThermo
    _WRITER_CONFIG = WriterConfig    
    _WRITER_CHECKPOINT = WriterCheckpoint

    def __init__(self, output_root, output_path, params, sim, swap_period, variables=['T']):
        self.output_root = output_root
        self.params = params
        self.variables = variables
        self.sim = sim
        self.steps_block = swap_period
        self.seed = 10
        self.nr = len(params)
        # The trajectory class is taken from the simulation backend
        # For the moment, we assume that at least the first simulation 
        # backend has a trajectory
        self.trajectory = self.sim[0].trajectory.__class__

        # Get physical replicas (systems) from simulation instances.
        # These are references: they'll follow the simulations 
        self.replica = [s.system for s in sim]

        # If output is just one directory, we listify it padding it with None
        if not isinstance(output_path, list):
            output_path = [output_path] + [None] * (self.nr-1)

        # Sanity check
        if not (self.nr == len(output_path) == len(sim)):
            raise ValueError('nr, params and sim must have the same len (%d, %d, %d)' % 
                             (self.nr, len(output_path), len(sim)))

        # Here it is good to call the base constructor since we know input sample
        # and output directory
        # TODO: potential bug here. If the initial system is different for each replica, the RMSD will be wrong
        # We should outsource rmsd or override.
        Simulation.__init__(self, self.replica[0], output_path)

        # Set random seed now
        random.seed(self.seed)

        self.replica_id = range(self.nr)
        # TODO: remember to checkpoint them
        self.accepted = [0.0 for i in range(self.nr)]
        self.attempts = [0.0 for i in range(self.nr)]
        # This is used to swap odd/even replicas in turn
        self.offset = 0

        # Distribute physical replicas in parallel.
        # Each process gets a bunch of replicas to evolve
        # my_replica contains their state id's.
        # We could as well distributes states, which would allow
        # for other optimizations.
        self._thermostat = None
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
        if not type(swap_period) is list:
            self.steps_block = [swap_period] * self.nr
        # If we provided 1 block period in list, upgrade it to nr length
        if len(self.steps_block) != self.nr:
            if len(self.steps_block) == 1:
                self.steps_block = swap_period * self.nr
            else:
                raise ValueError('number of blocks does not match n. replica')

        # We setup target steps here, then run_batch will use it
        # In principle, this could be controlled by the user from outside
        # and we could drop swap_period and steps_block completely
        # for s, steps in zip(self.sim, self.steps_block):
        #     s.setup(target_steps=steps)

        # Define output files
        mkdir(self.output_root + '/state')
        mkdir(self.output_root + '/replica')
        self.file_log = self.output_root + '/pt.log'
        # For each thermodynamic state, info on the replica which has it
        self.file_state_out = [self.output_root + '/state/%d.out' % i for i in range(self.nr)]
        self.file_state_xyz = [self.output_root + '/state/%d.xyz' % i for i in range(self.nr)]
        # For each physical replica, info on the state in which it is
        self.file_replica_out = [self.output_root + '/replica/%d.out' % i for i in range(self.nr)]
        # Make sure output directories exist
        for d in self.output_path:
            if d:
                mkdir(d)

    @property
    def rmsd(self):
        """In parallel tempering simulation we define rmsd as the minimum one"""
        return min([s.rmsd for s in self.sim])

    def clean_files(self):
        for f in \
                self.file_state_out + \
                self.file_state_xyz + \
                self.file_replica_out:
            rmf(f)

    def write_log(self):
        """ Dump state parameters """
        with open(self.file_log, 'w') as fh:
            for i, T in enumerate(self.params):
                fh.write('%d %g\n' % (i, T))

    def write_replica(self, msd, step):
        """ Dump output info on a physical replica """
        # Loop over replicas
        for i in range(self.nr):
            with open(self.file_replica_out[i], 'a') as fh:
                # In which state is physical replica i ?
                fh.write('%d %d %d %g\n' % (self.steps, step[i], self.state[i], msd[i]))

    def write_state(self, u, step):
        """ Dump output info on a thermodynamic state """
        # TODO: we could write state atomically, which would allow parallelization
        # Loop over states       
        logging.debug('rx step=%s replicas(state)=%s' % (step[0], self.replica_id))

        for i in range(self.nr):
            with open(self.file_state_out[i], 'a') as fh:
                # Which replica is in state i? What is its energy?
                fh.write('%d %d %d %g\n' % (self.steps, step[i], self.replica_id[i], u[self.replica_id[i]])) #, self.acceptance(i)))

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
        logging.debug('write checkpoint %d' % self.steps)
        # Note: offset and step are redundant, since they are global
        for i in self.my_replica:
            with open(self.file_replica_out[i] + '.chk', 'w') as fh:
                fh.write('%d\n' % self.state[i])
                fh.write('%d\n' % self.steps)
                fh.write('%d\n' % self.offset)
            # TODO: write_checkpoint is not part of the official simulation interface, should it?
            self.sim[i].write_checkpoint()

    def run_pre(self):
        Simulation.run_pre(self)

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
                    logging.debug('pt restarting replica %d at state %d from step %d' % 
                                  (i, istate, self.steps))
        # Restarting is handled by the simulation instance.
        # The only glitch for now is that the checkpoint file ends up
        # in the directory corresponding to the initial state of the replica
        for i in self.my_replica:
            self.sim[i].run_pre()
        
        # Log RX info
        logging.info('rx with %d GPUs (rank=%d)' % (size, rank))
        self.write_log()

        if not self.restart:
            self.clean_files()

    def run_until(self, n):
        # Evolve over my physical replicas.
        logging.info('run until %d' % n)
        for i in self.my_replica:
            logging.debug('evolve %d' % i)
            # This will evolve physical replica[i] for the number
            # of steps prescribed for its state 
            n = self.sim[i].steps + self.steps_block[self.state[i]]
            self.sim[i].run_until(n)
        # Attempt to exchange replicas.
        # When swapping we must exchange their thermostats.
        self.exchange(self.replica)
            
    def run_end(self):
        for i in self.my_replica:
            fout = self.output_root + '/state/%d.xyz' % self.state[i]
            # TODO: this is not good. None is not accepted by other formats. Why needed here?
            # Where do we use this?
#            with self.trajectory(fout, 'w') as t:
#                t.write_sample(self.replica[i], None)
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

    def _swap(self, system, state, nn_state, r_i, r_j):
        tmp = self.replica_id[state]
        self.replica_id[state] = self.replica_id[nn_state]
        self.replica_id[nn_state] = tmp

    def _process_with(self, state):
        return self.process_with_replica[self.replica_id[state]]

    def exchange(self, system):
        logging.debug("exchange rx")
        if self.variables == ['T']:
            if self._thermostat is None:
                self._thermostat = [s.thermostat for s in system]
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
        self._update_replicas()

        # Update internal state of replicas / simulations
        # When swapping thermostats we imply the internal state should be swapped.
        # This is not the case in RUMD for now (we simply swap T), so we simply 
        # communicate the state and use a static mapping between states and temperatures.
        # The internal state of the thermostat after each swap is reset to the initial one
        for i, s in enumerate(self.state):
            system[i].thermostat = self._thermostat[s]

    def _update_replicas(self):
        """Once states have changed, we must update replica id's across processes"""
        state_tmp = numpy.array(range(self.nr))
        state_new = numpy.array([self.state[r] for r in self.my_replica])
        if size > 1:
            comm.Allgather(state_new, state_tmp)
        else:
            state_tmp = state_new
        # After setting replica_id we can safely use self.state
        for i, s in enumerate(state_tmp):
            self.replica_id[s] = i

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
        T_i = self.params[my_state]
        T_j = self.params[nn_state]
        # Store current probability term
        log.debug("comm rank %d uj=%g uu=%g Ti=%g Tj=%g" % (rank,  u_j[0], u_i[0], T_i, T_j))
        self.__prob = math.exp(-(u_j[0]-u_i[0])*(1/T_i-1/T_j))
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

        # Final timer dump
#         if rank == 0:
#             fh = open(self.file_log, 'a')
#             fh.write("# Timings\n")
#             fh.write("\tNR \t NGPU \tTotal \tRun \tComm \tIO\n")
# #            sys.stdout.write("CPU \t%d \t%d \t%.2f \t%.2f \t%.2f \t%.2f\n" % (self.nr, size, t.cpu_time, t_run.cpu_time, t_com.cpu_time, t_io.cpu_time))
#             fh.write("WALL  \t%d \t%d \t%.2f \t%.2f \t%.2f \t%.2f\n" % (self.nr, size, t.wall_time, t_run.wall_time, t_com.wall_time, t_io.wall_time))
#             fh.close()


