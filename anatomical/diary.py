import sys
import logging

# Function to dynamically set up loggers in the bids_dir directory
def setup_logging(bids_dir):
    if not os.path.exists(bids_dir):
        os.makedirs(bids_dir)

    # Setup command logger
    command_logger = logging.getLogger("command_logger")
    command_logger.setLevel(logging.INFO)

    command_file = os.path.join(bids_dir, "command_diary.txt")
    command_handler = logging.FileHandler(command_file)
    command_formatter = logging.Formatter('%(asctime)s - COMMAND - %(message)s')
    command_handler.setFormatter(command_formatter)

    # Remove previous handlers if any, to avoid duplicating logs
    if command_logger.hasHandlers():
        command_logger.handlers.clear()

    command_logger.addHandler(command_handler)

    # Setup verbose logger
    verbose_logger = logging.getLogger("verbose_logger")
    verbose_logger.setLevel(logging.DEBUG)

    verbose_file = os.path.join(bids_dir, "verbose_output.txt")
    verbose_handler = logging.FileHandler(verbose_file)
    verbose_formatter = logging.Formatter('%(asctime)s - VERBOSE - %(message)s')
    verbose_handler.setFormatter(verbose_formatter)

    # Remove previous handlers if any, to avoid duplicating logs
    if verbose_logger.hasHandlers():
        verbose_logger.handlers.clear()

    verbose_logger.addHandler(verbose_handler)

    # Also print verbose output to the terminal
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(verbose_formatter)
    verbose_logger.addHandler(stream_handler)

    return command_logger, verbose_logger

# Trace function to log function calls and their arguments
def trace_calls(frame, event, arg, command_logger, verbose_logger):
    if event == 'call':
        func_name = frame.f_code.co_name
        if func_name == "preprocess_anat":  # Target the function you want to log
            args = frame.f_locals  # Access the arguments of the function
            command_logger.info(f"Called {func_name} with arguments: {args}")
            verbose_logger.debug(f"Running {func_name}...")
    elif event == 'return':
        func_name = frame.f_code.co_name
        if func_name == "preprocess_anat":
            verbose_logger.debug(f"{func_name} returned.")
    return lambda frame, event, arg: trace_calls(frame, event, arg, command_logger, verbose_logger)

# Enable tracing
def enable_tracing(command_logger, verbose_logger):
    sys.settrace(lambda frame, event, arg: trace_calls(frame, event, arg, command_logger, verbose_logger))