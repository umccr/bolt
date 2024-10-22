import logging
import sys
from pathlib import Path

def setup_logging(output_dir):
    log_file = Path(output_dir) / 'smvl_somatic.log'
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