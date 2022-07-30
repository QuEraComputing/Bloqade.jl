from .version import __version__

import juliacall
juliacall.Main.seval('using Bloqade')
juliacall.Main.seval('using PythonCall')

from .julia import *
