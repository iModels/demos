import os
import logging
import time

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

import mbuild as mb
import metamds as mds
import mdtraj as md


from simgen.project import Project

OFFLINE = True

def build_monolayer(chain_length, n_molecules, **kwargs):
    from mbuild.examples import AlkaneMonolayer
    pattern = mb.Random2DPattern(n_molecules)
    monolayer = AlkaneMonolayer(pattern, tile_x=1, tile_y=1, chain_length=chain_length)
    monolayer.name = 'alkane_n-{}_l-{}'.format(n_molecules, chain_length)
    mb.translate(monolayer, [0, 0, 2])
    box = monolayer.boundingbox
    monolayer.periodicity += np.array([0, 0, 5 * box.lengths[2]])

    return monolayer

def generate_code(**parameters):
    # import pdb; pdb.set_trace()

    if OFFLINE:
        # load simgen files from local folders
        res_dir = os.path.join(os.path.dirname(__file__), 'simgen_project')

        manifest = {
            'title': 'Monolayer example',
            'code_path': [os.path.join(res_dir, 'code')],
            'concept_path': [os.path.join(res_dir, 'concepts')],
            'template_path': [os.path.join(res_dir, 'templates')]
        }

        project = Project(manifest)
    else:
        # load sigmen files from GitHub
        project = Project(os.path.join(os.path.dirname(__file__), 'simgen_project', 'online_project.yaml'))

    return project.render_tasks('prg', output_dir='./', inject_dict=parameters)

if __name__ == '__main__':
    # # configure logging
    # logging.basicConfig(format='%(module)s:%(lineno)d %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=0)

    # Input parameters
    parameters = {'chain_length': 10,
                  'n_molecules': 100,
                  'forcefield': 'OPLS-aa',
                  }

    # Build the initial configuration
    compound = build_monolayer(**parameters)
    parameters['compound'] = compound
    parameters['system_name'] = compound.name

    #compound.visualize()

    # Initialize a simulation instance with a template and some metadata
    sim = mds.Simulation(name='monolayer', template=generate_code, input_dir='static_input_files', output_dir='simgen_output')

    # Parameterize our simulation template
    task = sim.parametrize(**parameters)

    # Run
    task.execute()
    # task.execute(hostname='rahman.vuse.vanderbilt.edu', username='ctk3b')
    # print(task.status())

    # time.sleep(10)
    # task.sync()

    # Analyze
    trajectories = task.get_output_files('trajectories')
    topologies = task.get_output_files('topologies')
    # Pick which one to select?

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