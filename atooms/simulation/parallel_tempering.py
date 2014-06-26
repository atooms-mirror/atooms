# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

import os
import sys
import math
import random
import copy
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

class ReplicaExchange(object):
    
    def __init__(self, params, variables=['T']):
        self.params = params
        self.variables = variables
        self.nr = len(params)
        self.replica = range(self.nr)
        self.accepted = [0.0 for i in range(self.nr)]
        self.attempts = [0.0 for i in range(self.nr)]
        # This is used to swap odd/even replicas in turn
        self.offset = 0

        # For parallel
        # TODO: use switch or something
        self._thermostat = None
        self.states = range(self.nr)
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

    # TODO: update state as soon as we swap
    @property
    def state(self):
        s = [0] * self.nr
        for i in range(self.nr):
            s[self.replica[i]] = i
        return s

    import copy

    def acceptance(self, state):
        if self.attempts[state] > 0:
            return self.accepted[state] / self.attempts[state]
        else:
            return 0.0

    def __repr__(self, state):
        return 'state %d <-> %d, replicas %d <-> %d, [%d/%d], prob %.4f]' % \
            (state, state+1, self.replica[state], self.replica[state+1], self.accepted[state], self.attempts[state], self.__prob)

    def _swap(self, system, state, nn_state, r_i, r_j):
        tmp = self.replica[state]
        self.replica[state] = self.replica[nn_state]
        self.replica[nn_state] = tmp

        # General method here?
        tmp = system[r_i].thermostat
        system[r_i].thermostat = system[r_j].thermostat
        system[r_j].thermostat = tmp

    def _process_with(self, state):
        return self.process_with_replica[self.replica[state]]

    def exchange(self, system):
        logging.debug("exchange rx")
        if self.variables == ['T']:
            if self._thermostat is None:
                self._thermostat = [s.thermostat for s in system]
            self.__exchange_T_parallel(system)
            # if size == 1:
            #     self.__exchange_T(system)
            # else:
            #     self.__exchange_T_parallel(system)
        else:
            raise ValueError(self.variables)

    def __exchange_T_parallel(self, system):

        # Loop over ALL states, find out which states are sending out messages
        barrier()
        sender = []
        focused = []
        for state in range(self.nr):
            if self.replica[state] in self.my_replica:
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

        # Update offset
        self.offset = (self.offset+1) % 2

        # When swapping thermostats we imply the internal state should be swapped.
        # This is not the case in RUMD for now (we simply swap T), so we simply 
        # communicate the state and use a static mapping between states and temperatures.
        # TODO: improve swapping of thermostat
        state_tmp = numpy.array(range(self.nr))
        state_new = numpy.array([self.state[r] for r in self.my_replica])
        if size > 1:
            comm.Allgather(state_new, state_tmp)
        else:
            state_tmp = state_new
        for i, s in enumerate(state_tmp):
            self.replica[s] = i
            system[i].thermostat = self._thermostat[s]
#        self.fh.write("comm rank %d all gather tmp %s \n" % (rank, state_tmp))

                    
    def __exchange_T_comm(self, my_state, nn_state, system):

        # TODO: add attempts
        exchange_attempt = False
        me = rank 
        # Compute potential energy
        r_i = self.replica[my_state]
        r_j = self.replica[nn_state] # This should be up to date
        u_i = numpy.array([system[r_i].potential_energy()])
        # Get process that has my nearest neighboring state
        nn = self._process_with(nn_state)
#        print "* comm rank", rank, my_state, "(replica %d)" % r_i,"->", nn_state, "nn", nn
#        self.fh.flush()
#        self.fh.write('nn %d rank %d\n' % (nn, rank))

        if nn == rank:
            if my_state > nn_state:
                return
#            self.fh.write('comm rank %d on same rank %d -> %d\n' % (rank, my_state, nn_state))
            u_j = numpy.array([system[r_j].potential_energy()])
            ran = numpy.array([random.random()])
        else:
            if my_state < nn_state:
#                self.fh.write('comm rank %d on /= rank %d -> %d\n' % (rank, my_state, nn_state))
                # I am on the left, I send first.
                u_j = numpy.array([0.0])
                comm.Send(u_i, nn, 10)
                comm.Recv(u_j, nn, 11)
                ran = numpy.array([random.random()])
                comm.Send(ran, nn, 12)
            elif my_state > nn_state:
#                self.fh.write('comm rank %d on /= rank %d <- %d\n' % (rank, my_state, nn_state))
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
#        self.fh.write("comm rank %d uj=%g uu=%g Ti=%g Tj=%g\n" % (rank,  u_j[0], u_i[0], T_i, T_j))
        self.__prob = math.exp(-(u_j[0]-u_i[0])*(1/T_i-1/T_j))
        # Test if we can swap states of replicas
#        self.fh.write("comm rank %d sawpping ran %g prob %g => %s\n" % (rank,  ran[0], self.__prob, ran[0]<self.__prob,))
        #logging.info("sawpping %g %g" % (self.__prob, ran[0]))
        if (ran[0] < self.__prob):
#            self.fh.write("comm rank %d sawpping successful \n" % (rank))
            self._swap(system, my_state, nn_state, r_i, r_j)
            # TODO: fix accepted counter

    def __exchange_T(self, system):
        """Simple non-parallel version"""
        # Sanity check
        if not len(system) == self.nr: raise ValueError

        # Loop over states defined by self.params with a jump of 2
        # The offset to swap only odd / even states at each attempt 
        for state in range(self.offset, self.nr-1, 2):
            self.attempts[state] += 1
            self.attempts[state+1] += 1
            # TODO: can we use thermostat temperatures directly here? Are params here just references to Ts? There is a sort duplication here no?
            T_i = self.params[state]
            T_j = self.params[state+1]
            # Index of physical replicas having states to swap
            r_i = self.replica[state]
            r_j = self.replica[state+1]
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


class WriterConfig(object):

    def __call__(self, e):
        logging.debug('writer config')
        # TODO: see everything belongs to rx except the trajectory, yes but it's the simulation responsibility to define writers etc
        for i in e.rx.my_replica:
            if e.rx_verbosity[e.rx.state[i]] > 0:
                e.trajectory[e.rx.state[i]].write_sample(e.replica[i], e.steps*e.steps_block[i], e.steps)

class WriterCheckpoint(object):

    def __call__(self, e):
        e.write_checkpoint()
        
class WriterThermo(object):

    def __call__(self, e):
        logging.debug('writer thermo')

        # Since we grab steps from simulations, we must gather them first
        # We could have each process write down its replicas and make it more efficient, see write_state()
        u = numpy.ndarray(e.rx.nr)
        u_l = numpy.ndarray(len(e.rx.my_replica))
        rmsd = numpy.ndarray(e.rx.nr)
        rmsd_l = numpy.ndarray(len(e.rx.my_replica))
        steps = numpy.ndarray(e.rx.nr, dtype=int)
        steps_l = numpy.ndarray(len(e.rx.my_replica), dtype=int)
        for i, ri in enumerate(e.rx.my_replica):
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

    def __init__(self, output, params, sim, swap_period, verbosity=None, variables=['T']):
        Simulation.__init__(self) # super does not work???
        
        self.output = output
        self.params = params
        self.variables = variables
        self.sim = sim
        self.steps_block = swap_period
        self.seed = 10
        self.trajectory = [s.trajectory for s in self.sim]
        random.seed(self.seed)

        # Get system objects from simulations
        # These are references: they'll follow the simulations 
        self.replica = [s.system for s in sim]
        # We also need a system, to adhere to the simulation interface
        # Essentially we need it for the RMSD at the moment, we should
        # think how to generalize TargetRMSD and rmsd(). Making the latter
        # a list might break other things in base class, it should be checked
        # Perhaps just the log, which should be encapsulated and would be
        # overridden here
        self.system = self.replica[0]

        self.rx = ReplicaExchange(params, variables)

        # Listify steps per block
        if not type(swap_period) is list:
            self.steps_block = [swap_period] * self.rx.nr

        # We setup target steps here, then run_batch will use it
        # In principle, this could be controlled by the user from outside
        # and we could drop swap_period and steps_block completely
        for s, steps in zip(self.sim, self.steps_block):
            s.setup(target_steps=steps)

        # Replica verbosity
        # TODO: handle list / scalar verbosities for RX
        self.rx_verbosity = [1] + [0] * (self.rx.nr-1)
        if verbosity:
            self.rx_verbosity = [verbosity] * self.rx.nr

        # Define output files
        mkdir(self.output + '/state')
        mkdir(self.output + '/replica')
        self.file_log = self.output + '/pt.log'
        self.file_state_out = [self.output + '/state/%d.out' % i for i in range(self.rx.nr)]
        self.file_replica_out = [self.output + '/replica/%d.out' % i for i in range(self.rx.nr)]

    def clean_files(self):
        for f in \
                self.file_state_out + \
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
        for i in range(self.rx.nr):
            with open(self.file_replica_out[i], 'a') as fh:
                # In which state is physical replica i ?
                fh.write('%d %d %d %g\n' % (step[i], self.rx.replica[self.rx.state[i]], self.rx.state[i], msd[i]))

    def write_state(self, u, step):
        """ Dump output info on a thermodynamic state """
        # TODO: we could write_state operate atomtically, which would allow parallelization
        # Loop over states       
        for i in range(self.rx.nr):
            logging.info('write_state %d' % step[i])
            with open(self.file_state_out[i], 'a') as fh:
                # Which replica is in state i? What is its energy?
                fh.write('%d %d %d %g\n' % (step[i], i, self.rx.replica[i], u[self.rx.replica[i]])) #, self.acceptance(i)))

    # def read_checkpoint(self):
    #     """ Checkpoint """
    #     for i in self.rx.my_replica:
    #         f = self.file_state_out[i] + '.chk'
    #         with open(f, 'r') as fh:
    #             self.rx.state[i] = fh.read()
    #         self.replica[i], steps = self.trajectory[i].read_checkpoint()

    def write_checkpoint(self):
        logging.debug('write checkpoint %d' % self.steps)
        for i in self.rx.my_replica:
            with open(self.file_state_out[i] + '.chk', 'w') as fh:
                fh.write('%d\n' % self.rx.state[i])
                fh.write('%d\n' % self.steps)
            self.sim[i].write_checkpoint()

    def run_pre(self):
        Simulation.run_pre(self)

        if self.restart:
            # TODO: steps should all be equal, we should check
            for i in self.rx.my_replica:
                # This is basically a read_checkpoint()
                f = self.file_state_out[i] + '.chk'
                if os.path.exists(f):
                    with open(f, 'r') as fh:
                        self.rx.state[i] = int(fh.readline())
                        self.steps = int(fh.readline())
                    logging.debug('restarting replica %d from step %d' % (i, self.steps))
                irx = self.rx.state[i]
                # Set restart in child simulations
                self.sim[irx].restart = True

        for i in self.rx.my_replica:
            irx = self.rx.state[i]
            self.sim[irx].run_pre()
        
        # Log RX info
        logging.info('rx with %d GPUs (rank=%d)' % (size, rank))
        self.write_log()

        if not self.restart:
            self.clean_files()

    def run_until(self, n):
        # Evolve all replicas. Loop is actually over states 
        # but we pick the current replica using replica list. 
        logging.info('n')
        for i in self.rx.my_replica:
            logging.debug('evolve %d' % i)
            # This will evolve state i, thus replica[i]
            n = (self.steps+1) * self.sim[i].target_steps
            self.sim[i].run_until(n)

        # Attempt to exchange replicas.
        # When swapping we must exchange their thermostats.
        self.rx.exchange(self.replica)
            
    def run_end(self):
        # Only close my files
        for i in self.rx.my_replica:
            self.trajectory[self.rx.state[i]].close()
        barrier()
        #t.stop()

        # Final timer dump
#         if rank == 0:
#             fh = open(self.rx.file_log, 'a')
#             fh.write("# Timings\n")
#             fh.write("\tNR \t NGPU \tTotal \tRun \tComm \tIO\n")
# #            sys.stdout.write("CPU \t%d \t%d \t%.2f \t%.2f \t%.2f \t%.2f\n" % (self.rx.nr, size, t.cpu_time, t_run.cpu_time, t_com.cpu_time, t_io.cpu_time))
#             fh.write("WALL  \t%d \t%d \t%.2f \t%.2f \t%.2f \t%.2f\n" % (self.rx.nr, size, t.wall_time, t_run.wall_time, t_com.wall_time, t_io.wall_time))
#             fh.close()


