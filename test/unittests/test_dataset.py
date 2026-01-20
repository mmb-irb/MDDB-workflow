import os
import pytest
import tempfile
import subprocess
import threading
import time
import multiprocessing
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
        ds = Dataset(dataset_path=str(db_path))
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
        ds.remove_entry(tmpdir / "wrong", verbose=True)
        assert ds.get_status(tmpdir / "wrong") is None  # 'wrong' project should be removed
        # Remove specific replicas from proj0
        ds.remove_entry(tmpdir / "proj0/replica0", verbose=True)
        assert ds.get_status(tmpdir / "proj0")['num_mds'] == 1
        ds.remove_entry(tmpdir / "proj1/replica0", verbose=True)
        ds.remove_entry(tmpdir / "proj1/replica1", verbose=True)
        assert ds.get_status(tmpdir / "proj1")['num_mds'] == 1  # One replica should remain

@pytest.mark.unit_int
def test_dataset_console_commands():
    """Test dataset console commands help output."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        db_path = tmpdir / "dataset.db"
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


# ============================================================================
# Lock file tests
# ============================================================================

@pytest.mark.unit_int
def test_database_lock_basic():
    """Test basic lock acquire and release."""
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = Path(tmpdir) / "test.db"
        db_path.touch()

        lock = DatabaseLock(db_path, timeout=5.0)
        lock_dir = db_path.with_stem('.lock_')

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


@pytest.mark.unit_int
def test_database_lock_reentrant():
    """Test that the same thread can acquire the lock multiple times."""
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = Path(tmpdir) / "test.db"
        db_path.touch()

        lock = DatabaseLock(db_path, timeout=5.0)
        lock_dir = db_path.with_stem('.lock_')

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


@pytest.mark.unit_int
def test_database_lock_context_manager():
    """Test lock context managers."""
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = Path(tmpdir) / "test.db"
        db_path.touch()

        lock = DatabaseLock(db_path, timeout=5.0)
        lock_dir = db_path.with_stem('.lock_')

        # Test write_lock context manager
        with lock.lock_file():
            assert lock_dir.exists()
            assert lock._lock_count == 1
        assert not lock_dir.exists()
        assert lock._lock_count == 0


@pytest.mark.unit_int
def test_database_lock_timeout():
    """Test that lock acquisition times out correctly."""
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = Path(tmpdir) / "test.db"
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


@pytest.mark.unit_int
def test_database_lock_force_release():
    """Test force release of stale locks."""
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = Path(tmpdir) / "test.db"
        db_path.touch()
        lock_dir = db_path.with_stem('.lock_')

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


@pytest.mark.unit_int
def test_database_lock_threading():
    """Test lock behavior with multiple threads."""
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = Path(tmpdir) / "test.db"
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


@pytest.mark.unit_int
def test_dataset_lock_integration():
    """Test that Dataset uses locks correctly for database operations."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        db_path = tmpdir / "dataset.db"

        # Create dataset
        ds = Dataset(dataset_path=str(db_path))

        # Verify lock directory is created during table creation
        lock_dir = db_path.with_stem('.lock_')
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


def _worker_update_status(db_path, uuid, project_uuid, worker_id, iterations, results_queue):
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


@pytest.mark.unit_int
def test_dataset_lock_multiprocess():
    """Test that locks work correctly across multiple processes."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        db_path = tmpdir / "dataset.db"

        # Create dataset and add a project
        ds = Dataset(dataset_path=str(db_path))
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
                target=_worker_update_status,
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
        lock_dir = db_path.with_stem('.lock_')
        assert not lock_dir.exists()

        # Verify the project still exists and has valid state
        ds = Dataset(dataset_path=str(db_path))
        status = ds.get_status(str(project_dir))
        assert status is not None
        assert status['state'] == State.RUNNING.value
        ds.close()
