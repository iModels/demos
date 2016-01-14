import itertools as it

import mbuild as mb
import metamds as mds
import mdtraj as md


def build_system(box, n_molecules, **kwargs):
    from mbuild.examples import Ethane
    ethane = Ethane()
    full_box = mb.fill_box(ethane, n_molecules, box)
    full_box.name = '{}_ethanes'.format(n_molecules)
    return full_box


def build_script(compound, forcefield, **kwargs):
    name = compound.name
    em = 'em.mdp'
    nvt = 'nvt.mdp'
    gro = '{name}.gro'.format(name=name)
    top = '{name}.top'.format(name=name)

    compound.save(top, forcefield=forcefield, overwrite=True)

    em_grompp = 'gmx grompp -f {mdp} -c {gro} -p {top} -o em.tpr'.format(
        mdp=em, gro=gro, top=top)
    em_mdrun = 'gmx mdrun -v -deffnm em'

    nvt_grompp = 'gmx grompp -f {mdp} -c em.gro -p {top} -o nvt.tpr'.format(
        mdp=nvt, top=top)
    nvt_mdrun = 'gmx mdrun -v -deffnm nvt'

    script = (em_grompp, em_mdrun, nvt_grompp, nvt_mdrun)
    return script


if __name__ == '__main__':
    sim = mds.Simulation(name='ethane', template=build_script, project_dir='ethane_box')

    # Input parameters
    parameters = {'n_molecules': 1,
                  'box': [3, 3, 3],
                  'forcefield': 'OPLS-aa'}

    # Build
    compound = build_system(**parameters)
    parameters['compound'] = compound

    task = sim.parametrize(**parameters)

    # Run
    job = task.execute(remote='', credentials='')

    exit()
    import ipdb; ipdb.set_trace()
    # Analyze
    job.status()

    output_traj = job.output['trajectories']
    # job.output['logs']
    traj = md.load(output_traj, top='')









