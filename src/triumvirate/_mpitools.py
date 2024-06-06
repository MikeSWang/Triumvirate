"""
MPI Tools (:mod:`~triumvirate._mpitools`)
==========================================================================

MPI-based parallelisation tools.

.. autosummary::
    allocate_tasks
    distribute_tasks

"""
import warnings

import numpy as np


def allocate_tasks(task_total, proc_total):
    """Allocate tasks to processes.

    This assigns ``tasks[i]`` tasks to the rank-``i`` process such that
    the total number of tasks, `task_total`, is shared amongst all
    `proc_total` processes.

    Parameters
    ----------
    task_total : int
        Total number of tasks.
    proc_total : int
        Total number of processes.

    Returns
    -------
    ntasks : list of int
        Number of tasks for each process.

    Raises
    ------
    TypeError
        When `task_total` or `proc_total` is not an integer.
    ValueError
        When `task_total` or `proc_total` is not positive.

    """
    if not isinstance(task_total, (int, np.integer)) \
            or not isinstance(proc_total, int):
        raise TypeError(
            "Both `task_total` and `proc_total` must be integers: "
            f"got {type(task_total)=} and {type(proc_total)=}"
        )
    if task_total <= 0 or proc_total <= 0:
        raise ValueError(
            "Both `task_total` and `proc_total` must be positive."
        )

    ntask_toassign, nproc_toassign, ntasks = task_total, proc_total, []
    while ntask_toassign > 0:
        ntask_assigned = ntask_toassign // nproc_toassign
        ntasks.append(ntask_assigned)
        ntask_toassign -= ntask_assigned
        nproc_toassign -= 1

    return ntasks


def distribute_tasks(ntasks=None, task_total=None, proc_total=None):
    """Distribute segments of tasks to each process.

    The number of tasks each process receives is determined by
    :func:`allocate_tasks` if `tasks` is `None`.  The rank-``i`` process
    receives ``ntasks[i]`` tasks in the range given by ``segments[i]``.

    Parameters
    ----------
    ntasks : list of int, optional
        Number of tasks each process receives.  If `None` (default),
        `task_total` and `proc_total` must be provided; otherwise,
        `task_total` and `proc_total` are both ignored.
    task_total : int, optional
        Total number of tasks (default is `None`).
    proc_total : int, optional
        Total number of processes (default is `None`).

    Returns
    -------
    segments : list of slice
        Index slices corresponding to the segment of tasks that each
        process should receive.

    """
    ntasks = ntasks if ntasks else allocate_tasks(task_total, proc_total)
    proc_total = proc_total if proc_total else len(ntasks)

    breakpoints = np.insert(np.cumsum(ntasks), 0, values=0)
    segments = [
        slice(breakpoints[rank], breakpoints[rank + 1])
        for rank in range(proc_total)
    ]

    return segments


def restore_warnings(captured_warnings, unique=True, show=True,
                     comm=None, comm_root=0):
    """Emit captured warnings.

    Parameters
    ----------
    captured_warnings : list of :class:`warnings.WarningMessage`
        List of recorded warnings as returned by
        ``warnings.catch_warnings(record=True)``.
    unique : bool, optional
        If `True` (default), only emit unique warnings; otherwise, emit
        all warnings.
    show : bool, optional
        If `True` (default), emit the warnings; otherwise, do not emit
        the warnings.
    comm : :class:`mpi4py.MPI.Comm`, optional
        MPI communicator (default is `None`).
    comm_root : int, optional
        Root process rank of the MPI communicator (default is 0).
        Ignored if `comm` is `None`.

    """
    if unique:
        restored_warnings = []
        seen_warnings = set()
        for record in captured_warnings:
            record_tuple = (
                record.message, record.category, record.filename, record.lineno
            )
            if record_tuple not in seen_warnings:
                restored_warnings.append(record)
                seen_warnings.add(record_tuple)
    else:
        restored_warnings = captured_warnings

    if show and (comm is None or comm.rank == comm_root):
        for record in restored_warnings:
            warnings.showwarning(
                record.message, record.category,
                record.filename, record.lineno,
                file=record.file, line=record.line
            )

    return restored_warnings
