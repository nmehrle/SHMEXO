#! /usr/bin/env python3
import os, shutil, argparse

parser = argparse.ArgumentParser()
parser.add_argument('--dir',
    default='../snap110f/',
    help='directory of the older version')
args = vars(parser.parse_args())

old_version = args['dir']
patch_file = old_version + 'patch_files'

with open(patch_file, 'r') as file:
  line = file.readline()
  # replacement files
  while len(line.strip()) > 0:
    if line[0] == '#':
      line = file.readline()
      continue
    new, _, old = line.split()
    os.system('mkdir -p %s' % os.path.dirname(new))
    #os.system('diff3 %s %s %s > %s.diff' % (old, old_version+old, old_version+new, new))
    print(old, old_version + new, new)
    os.system('diff %s %s > %s.diff' % (old_version + old, old_version + new, new))
    if not os.path.isfile(new):
      shutil.copy(old, new)
    line = file.readline()

  line = file.readline()
  # new files
  while len(line.strip()) > 0:
    if line[0] == '#':
      line = file.readline()
      continue
    fname = line.split()[0]
    os.system('mkdir -p %s' % os.path.dirname(fname))
    if not os.path.isfile(fname):
      shutil.copy(old_version+fname, fname+'.tmp')
    line = file.readline()
