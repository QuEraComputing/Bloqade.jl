import os
import toml
import click

root = os.path.dirname(__file__)
project_toml = os.path.join(root, 'pyproject.toml')

def read_project_toml():
    return toml.load(project_toml)

def update_project_toml(version):
    d = read_project_toml()
    d['tool']['poetry']['version'] = version
    with open(project_toml, 'w+') as f:
        toml.dump(d, f)

def write_version_file(version):
    with open(os.path.join(root, 'bloqade', 'version.py'), 'w+') as f:
        f.write('__version__ = "{}"\n'.format(version))

@click.command()
@click.option('--version', '-v', required=True, help='Version number to use.')
def bump_version(version):
    update_project_toml(version)
    write_version_file(version)
    return

if __name__ == '__main__':
    bump_version()
