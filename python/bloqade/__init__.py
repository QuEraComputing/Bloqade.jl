import juliacall
juliacall.Main.seval('using Bloqade')
juliacall.Main.seval('using PythonCall')


from .version import __version__
from .julia import *
