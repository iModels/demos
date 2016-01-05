import itertools as it

import mbuild as mb
import metamds as mds
import mdtraj as md


def build_system(n_particles, type1_fraction=0.5):
    compound = mb.Compound()
    return compound

if __name__ == '__main__':
    # YAML file defining sim tasks
    # hoomd input script
    # gromacs input script
    # lammps input script

    simulation = mds.Simulation(yml_file='')  # yml_file can be remote on github
    # simulation = mds.Simulation(github_repo='https...', credentials='')
        # happens internally in metamds if mbuild.compound is passed
        # compound.save('foo.top', forcefield=ff_file_handle)
        # compound.save('foo.hoomdxml', forcefield=ff_file_handle)
        # compound.save('foo.lmp', forcefield=ff_file_handle)

    ff_file_handle = ''

    temps = [1, 2, 3]
    fractions = [0.25, 0.5, 0.75]
    for temp, frac in it.product(temps, fractions):
        compound = build_system(100, frac)

        input_parameters = {'compound': compound,
                            'temperature': temp}
        task_name = 'lj_frac-{:.2f}_temp-{:d}'.format(frac, temp)
        task = simulation.prep(input_parameters, name=task_name)
        job = task.run(remote='', credentials='')

        job.status()

        output_traj = job.output['trajectories']
        # job.output['logs']
        traj = md.load(output_traj, top='')









