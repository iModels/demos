import time
import os

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

import mbuild as mb
import metamds as mds
import mdtraj as md


def build_monolayer(chain_length, n_molecules, **kwargs):
    from mbuild.examples import AlkaneMonolayer
    pattern = mb.Random2DPattern(n_molecules)
    monolayer = AlkaneMonolayer(pattern, tile_x=1, tile_y=1, chain_length=chain_length)
    monolayer.name = 'alkane_n-{}_l-{}'.format(n_molecules, chain_length)
    mb.translate(monolayer, [0, 0, 2])
    return monolayer


def create_run_script(build_func, forcefield, input_dir, **kwargs):
    compound = build_func(**kwargs)
    name = compound.name
    em = os.path.join(input_dir, 'em.mdp')
    nvt = os.path.join(input_dir, 'nvt.mdp')
    gro = '{name}.gro'.format(name=name)
    top = '{name}.top'.format(name=name)

    box = compound.boundingbox
    compound.periodicity += np.array([0, 0, 5 * box.lengths[2]])
    compound.save(top, forcefield=forcefield, overwrite=True)

    em_grompp = 'gmx grompp -f {mdp} -c {gro} -p {top} -o em.tpr'.format(mdp=em, gro=gro, top=top)
    em_mdrun = 'gmx mdrun -v -deffnm em -ntmpi 1'

    nvt_grompp = 'gmx grompp -f {mdp} -c em.gro -p {top} -o nvt.tpr'.format(mdp=nvt, top=top)
    nvt_mdrun = 'gmx mdrun -v -deffnm nvt -ntmpi 1'

    script = (em_grompp, em_mdrun, nvt_grompp, nvt_mdrun)
    return script

if __name__ == '__main__':
    # Input parameters
    parameters = {'chain_length': 10,
                  'n_molecules': 100,
                  'forcefield': 'OPLS-aa'}

    # Build the initial configuration
    compound = build_monolayer(**parameters)
    #compound.visualize()

    parameters['build_func'] = build_monolayer

    # Initialize a simulation instance with a template and some metadata
    sim = mds.Simulation(name='monolayer', template=create_run_script, output_dir='output')

    # Parameterize our simulation template
    task = sim.parametrize(**parameters)

    # Run
    # task.execute()
    # exit()
    task.execute(hostname='rahman.vuse.vanderbilt.edu', username='ctk3b')
    print(task.status())

    time.sleep(10)
    task.sync()

    # Analyze
    trajectories = task.get_output_files('trajectories')
    topologies = task.get_output_files('topologies')
    # Pick which one to select?
    import pdb; pdb.set_trace()

    trj_path = os.path.join(task.output_dir, 'nvt.xtc')
    top_path = os.path.join(task.output_dir, 'em.gro')
    traj = md.load(trj_path, top=top_path)
    print(traj)

    # RDF
    # pairs = traj.top.select_pairs('name C', 'name C')
    # r, g_r = md.compute_rdf(traj, pairs)
    # plt.plot(r, g_r)
    # plt.xlabel('r (nm)')
    # plt.ylabel('g(r)')
    # plt.show()
    #
    # s2 = md.compute_nematic_order(traj, 'residues')
    # plt.plot(traj.time, s2)
    # plt.xlabel('time (ps)')
    # plt.ylabel('S2')