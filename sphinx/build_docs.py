import os
import subprocess
import configparser
from pathlib import Path

MAINBRANCH = 'master'
BUILD_DIR = '../_docs'
PAGES_DIR = '../pages'
SPHINXSOURCE = './source'

def build_doc(version, tag):
    os.environ['current_version'] = version
    # checkout to the tagged commit
    subprocess.run(f'git checkout {tag}', shell=True)
    # Recover the latest conf.py and latest versions.yaml
    subprocess.run(f'git checkout {MAINBRANCH} -- {SPHINXSOURCE}/conf.py', shell=True)
    subprocess.run(f'git checkout {MAINBRANCH} -- {SPHINXSOURCE}/versions.ini', shell=True)
    # build the docs
    subprocess.run("make html", shell=True)

def move_dir(src, dst):
    print(f'    Creating: {Path(dst)}')
    Path(dst).mkdir(exist_ok=True)
    if not src.endswith('/'):
        src = f'{src}/'
    subprocess.run(f'mv {src}* {dst}', shell=True)

os.environ['build_all_docs'] = str(False)
os.environ['pages_root'] = 'https://maserasgroup-repo.github.io/pyssian' 

build_doc('latest', MAINBRANCH)
move_dir(f'{BUILD_DIR}/html/', f'{PAGES_DIR}/')

versions = configparser.ConfigParser(allow_no_value=True)
versions.read([f'{SPHINXSOURCE}/versions.ini'])

for version in versions.options('versions'):
    tag = versions['versions'][version]
    build_doc(version, tag)
    move_dir(f'{BUILD_DIR}/html/', f'{PAGES_DIR}/{version}/')