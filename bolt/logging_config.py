import logging
import sys
import pathlib
from datetime import datetime

class IgnoreTinfoFilter(logging.Filter):
    def filter(self, record):
        # Exclude messages that contain the unwanted text.
        if "no version information available" in record.getMessage():
            return False
        return True

def setup_logging(output_dir, script_name):
    # Create a timestamp for the log file
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_filename = f"{script_name}_{timestamp}.log"
    log_file = pathlib.Path(output_dir) / log_filename

    # Create individual handlers.
    console_handler = logging.StreamHandler(sys.stdout)
    file_handler = logging.FileHandler(log_file)

    # Instantiate and attach the filter to both handlers.
    tinfo_filter = IgnoreTinfoFilter()
    console_handler.addFilter(tinfo_filter)
    file_handler.addFilter(tinfo_filter)

    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[file_handler, console_handler]
    )
    logger = logging.getLogger(__name__)
    logger.info("Logging setup complete")
