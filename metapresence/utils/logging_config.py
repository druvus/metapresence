"""
Logging configuration for the Metapresence package.
"""

import logging
import sys


def setup_logging(verbose=False, quiet=False):
    """
    Set up logging configuration.
    
    Args:
        verbose: Enable verbose logging (DEBUG level).
        quiet: Minimize logging (only show WARN and above).
    """
    # Create logger
    logger = logging.getLogger('metapresence')
    
    # Set logging level based on verbosity
    if verbose:
        logger.setLevel(logging.DEBUG)
    elif quiet:
        logger.setLevel(logging.WARNING)
    else:
        logger.setLevel(logging.INFO)
    
    # Create console handler and set level
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.DEBUG)
    
    # Create formatter
    formatter = logging.Formatter('[%(levelname)s] %(message)s')
    
    # Add formatter to console handler
    console_handler.setFormatter(formatter)
    
    # Clear any existing handlers to avoid duplicates
    logger.handlers = []
    
    # Add console handler to logger
    logger.addHandler(console_handler)
    
    return logger


# Create a default logger
logger = logging.getLogger('metapresence')
