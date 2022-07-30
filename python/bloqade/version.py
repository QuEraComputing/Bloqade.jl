# workaround poetry version bump issue
# we have a long time loading Julia anyway
# so this is fine

import os
import toml
root = os.path.dirname(os.path.dirname(__file__))
d = toml.load(os.path.join(root, 'pyproject.toml'))
__version__ = d['tool']['poetry']['version']
