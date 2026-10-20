import os
import pytest
import tempfile
import subprocess
import threading
import time
import multiprocessing
import sqlite3
from mddb_workflow.core.dataset import Dataset, State, DatabaseLock
from pathlib import Path


@pytest.mark.unit_int
def test_add_remove_entries():
    """Test that adding projects works correctly."""
    # Create a temporary directory to act as the dataset root
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        # Create some fake project directories
        project_dirs = [tmpdir / f"proj{i}" for i in range(3)]
        # Add 2 replicas to the first project and 3 to the second
        replica_dirs = [project_dirs[0] / f"replica{i}" for i in range(2)]
        replica_dirs += [project_dirs[1] / f"replica{i}" for i in range(3)]
        project_dirs += [tmpdir / f"proj_glob{i}" for i in range(3)]
        project_dirs += [tmpdir / "wrong"]
        for d in project_dirs + replica_dirs:
            os.makedirs(d)
        db_path = tmpdir / "dataset.db"
        ds = Dataset(dataset_path=str(db_path), create=True)
        project_dirs = [str(d) for d in project_dirs]
        ds.add_entries(project_dirs,
                       ignore_dirs=str(tmpdir / "proj_glob2"),
                       md_dirs="replica*",
                       verbose=True)
        for d in project_dirs:
            status = ds.get_status(d)
            if "glob2" in d:
                assert status is None  # Ignored directory should not be added
            else:
                assert status is not None
                assert status['state'] == State.NEW.value
                assert status['message'] == 'No information recorded yet.'
            if "proj0" in d:
                assert status['num_mds'] == 2, f"Expected 2 MDs for {d}, got {status['num_mds']}"
            elif "proj1" in d:
                assert status['num_mds'] == 3
        # Now remove the 'wrong' project
        ds.remove_entry("wrong", verbose=True)
        assert ds.get_status(tmpdir / "wrong") is None  # 'wrong' project should be removed
        # Remove specific replicas from proj0
        assert ds.get_status(tmpdir / "proj0")['num_mds'] == 2
        ds.remove_entry("proj0/replica0", query_scope='mds', verbose=True)
        assert ds.get_status(tmpdir / "proj0")['num_mds'] == 1, "Expected 1 MD after removing replica0"
        ds.remove_entry("proj1/replica0", query_scope='mds', verbose=True)
        ds.remove_entry("proj1/replica1", query_scope='mds', verbose=True)
        assert ds.get_status(tmpdir / "proj1")['num_mds'] == 1  # One replica should remain


@pytest.mark.unit_int
def test_dataset_console_commands():
    """Test dataset console commands help output."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        # db_path = tmpdir / "dataset.db"

        # Test main dataset help
        result = subprocess.run(
            ['mwf', 'dataset', '--help'],
            capture_output=True,
            text=True
        )
        assert result.returncode == 0
        assert 'dataset' in result.stdout.lower()
        # Test add subcommand help
        result = subprocess.run(
            ['mwf', 'dataset', 'add', '--help'],
            capture_output=True,
            text=True
        )
        assert result.returncode == 0
        assert 'add' in result.stdout.lower()
        # Test show subcommand help
        result = subprocess.run(
            ['mwf', 'dataset', 'show', '--help'],
            capture_output=True,
            text=True
        )
        assert result.returncode == 0
        assert 'show' in result.stdout.lower()
        # Test inputs subcommand help
        result = subprocess.run(
            ['mwf', 'dataset', 'inputs', '--help'],
            capture_output=True,
            text=True
        )
        assert result.returncode == 0
        assert 'inputs' in result.stdout.lower()
        # Test run subcommand help
        result = subprocess.run(
            ['mwf', 'dataset', 'run', '--help'],
            capture_output=True,
            text=True
        )
        assert result.returncode == 0
        assert 'run' in result.stdout.lower()
        # Test scan subcommand help
        result = subprocess.run(
            ['mwf', 'dataset', 'scan', '--help'],
            capture_output=True,
            text=True
        )
        assert result.returncode == 0
        assert 'scan' in result.stdout.lower()
        # Test watch subcommand help
        result = subprocess.run(
            ['mwf', 'dataset', 'watch', '--help'],
            capture_output=True,
            text=True
        )
        assert result.returncode == 0
        assert 'watch' in result.stdout.lower()


@pytest.mark.unit_int
class TestDatabaseLock:
    """Test the DatabaseLock class."""

    db_name = "test_lock.db"

    def _paths(self, tmpdir):
        db_path = Path(tmpdir) / self.db_name
        lock_dir = DatabaseLock._path_to_lock_dir(db_path)
        return db_path, lock_dir

    def test_database_lock_basic(self):
        """Test basic lock acquire and release."""
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path, lock_dir = self._paths(tmpdir)
            db_path.touch()

            lock = DatabaseLock(db_path, timeout=5.0)
            # Lock should not exist initially
            assert not lock_dir.exists()
            assert not lock.is_locked()

            # Acquire lock
            lock.acquire()
            assert lock_dir.exists()
            assert lock_dir.is_dir()
            assert lock.is_locked()

            # Release lock
            lock.release()
            assert not lock_dir.exists()
            assert not lock.is_locked()

    def test_database_lock_reentrant(self):
        """Test that the same thread can acquire the lock multiple times."""
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path, lock_dir = self._paths(tmpdir)
            db_path.touch()

            lock = DatabaseLock(db_path, timeout=5.0)
            # Acquire lock multiple times (reentrant)
            lock.acquire()
            lock.acquire()
            lock.acquire()
            assert lock._lock_count == 3

            # Release once - lock should still be held
            lock.release()
            assert lock._lock_count == 2
            assert lock_dir.exists()

            # Release again
            lock.release()
            assert lock._lock_count == 1
            assert lock_dir.exists()

            # Final release - lock should be gone
            lock.release()
            assert lock._lock_count == 0
            assert not lock_dir.exists()

    def test_database_lock_context_manager(self):
        """Test lock context managers."""
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path, lock_dir = self._paths(tmpdir)
            db_path.touch()

            lock = DatabaseLock(db_path, timeout=5.0)
            # Test write_lock context manager
            with lock.lock_file():
                assert lock_dir.exists()
                assert lock._lock_count == 1
            assert not lock_dir.exists()
            assert lock._lock_count == 0

    def test_database_lock_timeout(self):
        """Test that lock acquisition times out correctly."""
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path, _ = self._paths(tmpdir)
            db_path.touch()

            lock1 = DatabaseLock(db_path, timeout=5.0)
            lock2 = DatabaseLock(db_path, timeout=0.3, retry_interval=0.05)

            # First lock acquires
            lock1.acquire()

            # Second lock should timeout
            start = time.time()
            with pytest.raises(TimeoutError) as exc_info:
                lock2.acquire()
            elapsed = time.time() - start

            assert "Could not acquire" in str(exc_info.value)
            assert elapsed >= 0.3
            assert elapsed < 1.0  # Should not take too long

            lock1.release()

    def test_database_lock_force_release(self):
        """Test force release of stale locks."""
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path, lock_dir = self._paths(tmpdir)
            db_path.touch()

            # Simulate a stale lock (directory exists but no process holds it)
            os.mkdir(lock_dir)
            assert lock_dir.exists()

            lock = DatabaseLock(db_path, timeout=0.1)

            # Normal acquire should fail/timeout
            with pytest.raises(TimeoutError):
                lock.acquire()

            # Force release should work
            lock.force_release()
            assert not lock_dir.exists()

            # Now acquire should work
            lock.acquire()
            lock.release()

    def test_database_lock_threading(self):
        """Test lock behavior with multiple threads."""
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path, _ = self._paths(tmpdir)
            db_path.touch()

            counter = {'value': 0}
            errors = []

            def increment_with_lock():
                try:
                    lock = DatabaseLock(db_path, timeout=10.0, retry_interval=0.01)
                    for _ in range(10):
                        with lock.lock_file():
                            current = counter['value']
                            time.sleep(0.001)  # Small delay to increase contention
                            counter['value'] = current + 1
                except Exception as e:
                    errors.append(str(e))

            # Start multiple threads
            threads = [threading.Thread(target=increment_with_lock) for _ in range(5)]
            for t in threads:
                t.start()
            for t in threads:
                t.join()

            # Check no errors occurred
            assert len(errors) == 0, f"Errors occurred: {errors}"
            # Check counter is correct (5 threads * 10 increments = 50)
            assert counter['value'] == 50

    def test_dataset_lock_integration(self):
        """Test that Dataset uses locks correctly for database operations."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            db_path, lock_dir = self._paths(tmpdir)

            # Create dataset
            ds = Dataset(dataset_path=str(db_path), create=True)

            # Lock should be released after __init__
            assert not lock_dir.exists()

            # Create a project directory
            project_dir = tmpdir / "proj1"
            project_dir.mkdir()

            # Add project (should use write lock internally)
            ds.add_project(str(project_dir), make_uuid=True)

            # Lock should be released after operation
            assert not lock_dir.exists()

            # Query status (should use read lock internally)
            status = ds.get_status(str(project_dir))
            assert status is not None
            assert status['state'] == State.NEW.value

            # Lock should be released after operation
            assert not lock_dir.exists()

            ds.close()

    def _worker_update_status(self, db_path, uuid, project_uuid, worker_id, iterations, results_queue):
        """Worker function for multiprocessing test."""
        try:
            ds = Dataset(dataset_path=str(db_path), lock_timeout=30.0)
            for i in range(iterations):
                ds.update_status(
                    uuid=uuid,
                    state=State.RUNNING,
                    message=f"Worker {worker_id} iteration {i}",
                    project_uuid=project_uuid if project_uuid else None,
                )
                time.sleep(0.01)  # Small delay
            ds.close()
            results_queue.put(('success', worker_id))
        except Exception as e:
            results_queue.put(('error', worker_id, str(e)))

    def test_dataset_lock_multiprocess(self):
        """Test that locks work correctly across multiple processes."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            db_path, lock_dir = self._paths(tmpdir)

            # Create dataset and add a project
            ds = Dataset(dataset_path=str(db_path), create=True)
            project_dir = tmpdir / "proj1"
            project_dir.mkdir()
            ds.add_project(str(project_dir), make_uuid=True)

            # Get the UUID for the project
            status = ds.get_status(str(project_dir))
            uuid = status['uuid']
            ds.close()

            # Start multiple processes that update the same project
            num_workers = 4
            iterations = 5
            results_queue = multiprocessing.Queue()

            processes = []
            for i in range(num_workers):
                p = multiprocessing.Process(
                    target=self._worker_update_status,
                    args=(db_path, uuid, None, i, iterations, results_queue)
                )
                processes.append(p)
                p.start()

            # Wait for all processes to complete
            for p in processes:
                p.join(timeout=60)

            # Collect results
            results = []
            while not results_queue.empty():
                results.append(results_queue.get())

            # Check all workers succeeded
            errors = [r for r in results if r[0] == 'error']
            assert len(errors) == 0, f"Worker errors: {errors}"
            assert len(results) == num_workers

            # Verify the lock directory is cleaned up
            assert not lock_dir.exists()

            # Verify the project still exists and has valid state
            ds = Dataset(dataset_path=str(db_path))
            status = ds.get_status(str(project_dir))
            assert status is not None
            assert status['state'] == State.ERROR.value
            assert 'no longer running' in status['message'].lower()
            ds.close()


def _set_running_status_from_child(db_path: str, project_dir: str, ready_queue=None, hold_event=None):
    """Set a project status to RUNNING from a child process."""
    ds = Dataset(dataset_path=str(db_path))
    status = ds.get_status(str(project_dir))

    ds.update_status(
        uuid=status['uuid'],
        state=State.RUNNING,
        message='Running workflow...',
    )

    if ready_queue is not None:
        ready_queue.put('ready')
    if hold_event is not None:
        hold_event.wait(timeout=30)
    ds.close()


@pytest.mark.unit_int
def test_running_project_lifecycle_reconciles_to_error_when_process_dies():
    """Project stays RUNNING while child is alive, then reconciles to ERROR after exit."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        db_path = tmpdir / "dataset.db"
        project_dir = tmpdir / "proj1"
        project_dir.mkdir()

        ds = Dataset(dataset_path=str(db_path), create=True)
        ds.add_project(str(project_dir), make_uuid=True)
        ready_queue = multiprocessing.Queue()
        hold_event = multiprocessing.Event()
        child = multiprocessing.Process(
            target=_set_running_status_from_child,
            args=(str(db_path), str(project_dir), ready_queue, hold_event),
        )
        child.start()
        assert ready_queue.get(timeout=10) == 'ready', "Child process did not signal ready in time"

        _ = ds.get_dataframe(query_scope='project')
        updated = ds.get_status(str(project_dir))
        assert updated is not None
        assert updated['state'] == State.RUNNING.value

        hold_event.set()  # Allow child to exit
        child.join(timeout=10)  # Wait for child to exit
        assert child.exitcode == 0

        # Once process exits, next refresh should reconcile RUNNING -> ERROR.
        _ = ds.get_dataframe(query_scope='project')
        updated = ds.get_status(str(project_dir))
        assert updated is not None
        assert updated['state'] == State.ERROR.value
        assert 'no longer running' in updated['message'].lower()
        ds.close()


def _open_dataset_concurrently(db_path: str, start_event, results_queue, iterations: int = 5):
    """Open the dataset several times and read from it, reporting any error."""
    try:
        start_event.wait(timeout=30)
        for _ in range(iterations):
            ds = Dataset(dataset_path=db_path, create=True)
            ds.query_table('projects')
            ds.close()
        results_queue.put(('success', os.getpid()))
    except Exception as e:
        results_queue.put(('error', os.getpid(), repr(e)))


def _run_concurrent_opens(db_path: Path, num_workers: int = 8) -> list:
    start_event = multiprocessing.Event()
    results_queue = multiprocessing.Queue()
    processes = [
        multiprocessing.Process(target=_open_dataset_concurrently, args=(str(db_path), start_event, results_queue))
        for _ in range(num_workers)
    ]
    for p in processes:
        p.start()
    start_event.set()  # Release all workers at the same time
    results = [results_queue.get(timeout=120) for _ in processes]
    for p in processes:
        p.join(timeout=30)
    return results


@pytest.mark.unit_int
def test_dataset_concurrent_open():
    """Several processes opening the same (new and existing) dataset at once all succeed."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        db_path = tmpdir / "dataset.db"

        # All workers race to create the schema of a brand new dataset
        results = _run_concurrent_opens(db_path)
        errors = [r for r in results if r[0] == 'error']
        assert not errors, f"Worker errors: {errors}"

        ds = Dataset(dataset_path=str(db_path))
        project_dir = tmpdir / "proj1"
        project_dir.mkdir()
        ds.add_project(str(project_dir), make_uuid=True)
        ds.close()

        # Opening an up to date dataset must not write to it
        before = db_path.read_bytes()
        mtime_before = db_path.stat().st_mtime_ns
        results = _run_concurrent_opens(db_path)
        errors = [r for r in results if r[0] == 'error']
        assert not errors, f"Worker errors: {errors}"
        assert db_path.read_bytes() == before
        assert db_path.stat().st_mtime_ns == mtime_before


@pytest.mark.unit_int
def test_dataset_open_retries_transient_errors(monkeypatch):
    """Transient sqlite3 errors while opening are retried; permanent ones are raised at once."""
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = Path(tmpdir) / "dataset.db"
        real_connect = sqlite3.connect
        failures = [sqlite3.OperationalError('unable to open database file'),
                    sqlite3.DatabaseError('database disk image is malformed')]
        calls = []

        def flaky_connect(*args, **kwargs):
            calls.append(None)
            if failures:
                raise failures.pop(0)
            return real_connect(*args, **kwargs)

        monkeypatch.setattr(Dataset, 'open_retry_base_delay', 0.01)
        monkeypatch.setattr(sqlite3, 'connect', flaky_connect)
        ds = Dataset(dataset_path=str(db_path), create=True)
        assert len(calls) == 3
        ds.close()

        def broken_connect(*args, **kwargs):
            calls.append(None)
            raise sqlite3.OperationalError('no such table: foo')

        calls.clear()
        monkeypatch.setattr(sqlite3, 'connect', broken_connect)
        with pytest.raises(sqlite3.OperationalError):
            Dataset(dataset_path=str(db_path))
        assert len(calls) == 1


@pytest.mark.unit_int
def test_lock_dir_next_to_writable_dataset():
    """The lock lives next to the database so every node of a shared filesystem sees it."""
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = Path(tmpdir) / "dataset.db"
        lock_dir = DatabaseLock._path_to_lock_dir(db_path)
        assert lock_dir == db_path.resolve().with_name('.lock_dataset.db')

        ds = Dataset(dataset_path=str(db_path), create=True)
        assert ds._lock.lock_dir == lock_dir
        with ds.locked_storage_file:
            assert lock_dir.is_dir()
        assert not lock_dir.exists()
        ds.close()


@pytest.mark.unit_int
@pytest.mark.skipif(hasattr(os, 'geteuid') and os.geteuid() == 0, reason="root can write to read-only directories")
def test_lock_dir_falls_back_to_tmp_for_read_only_dataset():
    """A dataset in a read-only directory (e.g. a coworker's) uses a lock in the temp dir and can still be read."""
    with tempfile.TemporaryDirectory() as tmpdir:
        data_dir = Path(tmpdir) / "coworker"
        data_dir.mkdir()
        db_path = data_dir / "dataset.db"
        project_dir = data_dir / "proj1"
        project_dir.mkdir()
        ds = Dataset(dataset_path=str(db_path), create=True)
        ds.add_project(str(project_dir), make_uuid=True)
        ds.close()

        data_dir.chmod(0o555)
        db_path.chmod(0o444)
        try:
            lock_dir = DatabaseLock._path_to_lock_dir(db_path)
            assert lock_dir.parent == Path(tempfile.gettempdir())
            assert lock_dir.name.startswith('.lock_dataset.db_')

            ds = Dataset(dataset_path=str(db_path))
            assert ds._lock.lock_dir == lock_dir
            assert len(ds.query_table('projects')) == 1
            assert not lock_dir.exists()
            ds.close()
        finally:
            db_path.chmod(0o644)
            data_dir.chmod(0o755)
