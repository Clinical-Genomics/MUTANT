import logging
import json
import os

__version__ = "0.1"

# Keep track of microSALT installation
wd = os.path.dirname(os.path.realpath(__file__))

# Load configuration
preset_config = ""
logger = ""
default = os.path.join(os.environ["HOME"], ".MUTANT/config.json")

if "MUTANT_CONFIG" in os.environ:
    try:
        envvar = os.environ["MUTANT_CONFIG"]
        with open(envvar, "r") as conf:
            preset_config = json.load(conf)
    except Exception as e:
        print("Config error: {}".format(str(e)))
        pass
elif os.path.exists(default):
    try:
        with open(os.path.abspath(default), "r") as conf:
            preset_config = json.load(conf)
    except Exception as e:
        print("Config error: {}".format(str(e)))
        pass

# Config dependent section:
if preset_config != "":
    try:

        # Initialize logger
        logger = logging.getLogger("main_logger")
        logger.setLevel(logging.INFO)
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        ch.setFormatter(logging.Formatter("%(levelname)s - %(message)s"))
        logger.addHandler(ch)

        fh = logging.FileHandler(
            os.path.expanduser(preset_config["folders"]["log_file"])
        )
        fh.setFormatter(
            logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        )
        logger.addHandler(fh)

    except Exception as e:
        print("Config error: {}".format(str(e)))
        pass