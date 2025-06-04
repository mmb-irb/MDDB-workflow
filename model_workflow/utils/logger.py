import logging
import sys
import io

class StreamToLogger:
    """
    Fake file-like stream object that redirects writes to a logger instance.
    """
    def __init__(self, logger, log_level=logging.INFO):
        self.logger = logger
        self.log_level = log_level
        self.linebuf = ''

    def write(self, buf):
        for line in buf.rstrip().splitlines():
            if line: # Only log non-empty lines
                self.logger.log(self.log_level, line.strip())

    def flush(self):
        # This is important for some scenarios, but for simple print() it might not be strictly necessary
        pass

def init_logger(log_path: str = 'output.log', slurm: bool = False):
    if slurm:
        # Make the system output stream to not be buffered
        # This is useful to make prints work on time in Slurm
        # Otherwise, output logs are written after the script has fully run
        # Note that this fix affects all modules and built-ins
        unbuffered = io.TextIOWrapper(open(sys.stdout.fileno(), 'wb', 0), write_through=True)
        sys.stdout = unbuffered

    # Configure a root logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO) # Or DEBUG

    # Create a file handler to save logs to a file
    file_handler = logging.FileHandler(log_path, mode='a')  # 'w' to overwrite, 'a' to append
    file_handler.setLevel(logging.INFO) # Or DEBUG
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # Log to console (pytest usually captures this anyway)
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging.INFO)
    logger.addHandler(stream_handler)
    
    sys.stdout = StreamToLogger(logger, logging.INFO)  # Redirect print to logger
    sys.stderr = StreamToLogger(logger, logging.ERROR)  # Redirect errors to logger

    return logger


# Test the logger when this file is run directly
if __name__ == "__main__":
    init_logger()
    print("This is a test print message that should be captured by the logger")
    logging.info("This is a direct logger.info message")
    logging.warning("This is a warning message")
    logging.error("This is an error message")
    
    try:
        raise ValueError("This is a test exception")
    except Exception as e:
        logging.exception("Caught an exception: %s", str(e))
    
    print("Check the output.log file for the logged messages")
