from pathlib import Path


def inputs_generator(project_dir: str):
    """Discover the input files of a GROMACS nanobody project from its folder structure.

    Layout (one project = one nanobody run):
        <project>/<name>_solv_ions.gro   -> solvated structure (matches the trajectory atoms)
        <project>/r1|r2|r3/md.xtc        -> one production trajectory per replica
        <project>/r1|r2|r3/md.tpr        -> per-replica topology (self-contained)
    """
    project_path = Path(project_dir)

    # Replicas: any subdir named r<number> that actually holds a production md.xtc, sorted r1, r2, r3...
    replica_dirs = sorted(
        (d for d in project_path.glob('r[0-9]*') if d.is_dir() and (d / 'md.xtc').exists()),
        key=lambda d: int(d.name[1:]),
    )
    if not replica_dirs:
        raise ValueError(f"No replica directories with md.xtc found in {project_path}")

    mds = [
        {
            'mdir': d.name,                                   # real folder name, e.g. 'r1'
            'traj': (d / 'md.xtc').relative_to(project_path).as_posix(),
        }
        for d in replica_dirs
    ]

    # Topology: a .tpr is self-contained (topology + coordinates); use the first replica's production
    # topology as the representative project topology (all replicas share the same molecular topology).
    top = (replica_dirs[0] / 'md.tpr').relative_to(project_path).as_posix()

    # Structure: the solvated + ionised structure, which matches the solvated trajectory atom set.
    solv = sorted(project_path.glob('*_solv_ions.gro'))
    struct = solv[0].relative_to(project_path).as_posix() if solv else None

    return {'mds': mds, 'top': top, 'struct': struct}


if __name__ == '__main__':
    import sys, json
    print(json.dumps(inputs_generator(sys.argv[1]), indent=2))
