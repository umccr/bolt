import logging
import sys
import pathlib


def setup_logging(output_dir, log_filename):
    log_file = pathlib.Path(output_dir) / log_filename
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    logger = logging.getLogger(__name__)
    logger.info("Logging setup complete")