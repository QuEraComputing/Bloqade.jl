# workaround poetry version bump issue
# we have a long time loading Julia anyway
# so this is fine

from asyncore import write
import os
import toml

def build(setup_kwargs):
    root = os.path.dirname(__file__)
    d = toml.load(os.path.join(root, 'pyproject.toml'))
    __version__ = d['tool']['poetry']['version']
    with open(os.path.join(root, 'bloqade', 'version.py'), 'w+') as f:
        f.write('__version__ = "{}"\n'.format(__version__))
