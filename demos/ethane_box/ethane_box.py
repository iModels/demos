import os
import time

import matplotlib.pyplot as plt
import seaborn as sns

import mbuild as mb
import metamds as mds
import mdtraj as md


def build_ethane_box(box, n_molecules, **kwargs):
    from mbuild.examples import Ethane
    ethane = Ethane()
    full_box = mb.fill_box(ethane, n_molecules, box)
    full_box.name = '{}_ethanes'.format(n_molecules)
    return full_box


def create_run_script(compound, forcefield, input_dir, **kwargs):
    name = compound.name
    em = os.path.join(input_dir, 'em.mdp')
    nvt = os.path.join(input_dir, 'nvt.mdp')
    gro = os.path.join(input_dir, '{name}.gro'.format(name=name))
    top = os.path.join(input_dir, '{name}.top'.format(name=name))

    compound.save(top, forcefield=forcefield, overwrite=True)

    em_grompp = 'gmx grompp -f {mdp} -c {gro} -p {top} -o em.tpr'.format(mdp=em, gro=gro, top=top)
    em_mdrun = 'gmx mdrun -v -deffnm em'

    nvt_grompp = 'gmx grompp -f {mdp} -c em.gro -p {top} -o nvt.tpr'.format(mdp=nvt, top=top)
    nvt_mdrun = 'gmx mdrun -v -deffnm nvt'

    script = (em_grompp, em_mdrun, nvt_grompp, nvt_mdrun)
    return script

if __name__ == '__main__':
    # Input parameters
    parameters = {'n_molecules': 200,
                  'box': [3, 3, 3],
                  'forcefield': 'OPLS-aa'}

    # Build the initial configuration
    compound = build_ethane_box(**parameters)
    #compound.visualize()

    parameters['compound'] = compound

    # Initialize a simulation instance with a template and some metadata
    sim = mds.Simulation(name='ethane', template=create_run_script, output_dir='output')

    # Parameterize our simulation template
    task = sim.parametrize(**parameters)

    # Run
    #task.execute()
    task.execute(hostname='rahman.vuse.vanderbilt.edu', username='ctk3b')
    print(task.status())

    time.sleep(10)
    task.sync()

    # Analyze
    trajectories = task.get_output_files('trajectories')
    topologies = task.get_output_files('topologies')
    # Pick which one to select?

    trj_path = os.path.join(task.output_dir, 'nvt.xtc')
    top_path = os.path.join(task.output_dir, 'em.gro')
    traj = md.load(trj_path, top=top_path)
    print(traj)
    import pdb; pdb.set_trace()

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
