import logging

version = '3.4.1'

# Initialize log
log = logging.getLogger("main_log")
log.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(logging.Formatter("[%(asctime)s] - %(levelname)s: %(message)s"))
log.addHandler(ch)
