
from .compute_tad_representation import compute_TAD_signature, to_log_odds
from .tadmap_base_config import tadmap_loglevel

__all__ = ['compute_TAD_signature', 'to_log_odds', 'tadmap_loglevel']


import pkgutil, pathlib, importlib

# from pkgutil import iter_modules
# from pathlib import Path
# from importlib import import_module

# https://julienharbulot.com/python-dynamical-import.html
# iterate through the modules in the current package
#
package_dir = pathlib.Path(__file__).resolve().parent
for (_, module_name, _) in pkgutil.iter_modules([package_dir]):
    if 'data' in module_name:
        module = importlib.import_module(f"{__name__}.{module_name}")
