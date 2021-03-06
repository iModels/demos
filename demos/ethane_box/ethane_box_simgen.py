import logging
import os
import time

import matplotlib.pyplot as plt
import seaborn as sns

import mbuild as mb
import metamds as mds
import mdtraj as md

from simgen.project import Project

OFFLINE = True

def build_ethane_box(box, n_molecules, **kwargs):
    from mbuild.examples import Ethane
    ethane = Ethane()
    full_box = mb.fill_box(ethane, n_molecules, box)
    full_box.name = '{}_ethanes'.format(n_molecules)
    return full_box

def generate_code(**parameters):
    # import pdb; pdb.set_trace()

    if OFFLINE:
        # load simgen files from local folders
        res_dir = os.path.join(os.path.dirname(__file__), 'binary_lj_sim')

        manifest = {
            'title': 'Binary LJ Simulation Test Project with mBuild',
            'code_path': [os.path.join(res_dir, 'code')],
            'concept_path': [os.path.join(res_dir, 'concepts')],
            'template_path': [os.path.join(res_dir, 'templates')]
        }

        project = Project(manifest)
    else:
        # load sigmen files from GitHub
        project = Project(os.path.join(os.path.dirname(__file__), 'binary_lj_sim', 'online_project.yaml'))

    return project.render_tasks('prg', output_dir='./', inject_dict=parameters)

if __name__ == '__main__':
    # # configure logging
    # logging.basicConfig(format='%(module)s:%(lineno)d %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=0)

    # Build the initial configuration
    ethane_box = build_ethane_box(n_molecules=200, box=[3, 3, 3])
    # ethane_box.visualize()


    # Input parameters
    parameters = {'compound': ethane_box,
                  'forcefield': 'OPLS-aa',
                  'system_name': 'ethane_box'}

    # Initialize a simulation instance with a template and some metadata
    sim = mds.Simulation(name='ethane', template=generate_code, input_dir='static_input_files', output_dir='simgen_output')

    # Parameterize our simulation template
    task = sim.parametrize(**parameters)

    # import pdb; pdb.set_trace()

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
    # import pdb; pdb.set_trace()

    trj_path = os.path.join(task.output_dir, 'nvt.xtc')
    top_path = os.path.join(task.output_dir, 'em.gro')
    traj = md.load(trj_path, top=top_path)
    print(traj)
    # import pdb; pdb.set_trace()

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
