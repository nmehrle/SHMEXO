#! /usr/bin/env python3
import fnmatch
import os, re
import subprocess

# recursively obtain all source files in a directory
def SourceFiles(directory, pattern):
  matches = []
  for root, dirnames, filenames in os.walk(directory):
    for filename in fnmatch.filter(filenames, pattern):
      if not 'z.junk' in root:
        matches.append(os.path.join(root, filename))
  return matches

# generate patch file
def PatchStructure(patch):
  # athena source files
  files = SourceFiles('src', '*.?pp')
  athena = [f[4:] for f in files]

  replacement, addition = [], []
  files = SourceFiles(patch, '*.?pp')
  for fname in files:
    source = fname[len(patch)+1:]
    if source in athena:
      replacement.append(source)
    else:
      addition.append(source)
  return replacement, addition

# remove all softlinks
def CleanSoftlinks(directory):
  os.system('find -L %s -xtype l -delete' % directory)
  files = SourceFiles('src', '*.?pp.old')
  for f in files:
    if not os.path.isfile(f[:-4]):
      os.rename(f, f[:-4])
    else:
      raise SystemExit('File %s exists' % f[:-4])

# write patch file
def WritePatchFile(fname, patches):
  CleanSoftlinks('src')
  files = []
  if os.path.isfile(fname):
    os.remove(fname)
  for patch in patches[::-1]:
    replacement, addition = PatchStructure(patch)
    with open(fname, 'w+') as file:
      file.write('# Replacement files:\n')
      for f in replacement:
        if f not in files:
          file.write('%s/%s -> src/%s\n' % (patch, f, f))
          files.append(f)
      file.write('\n# Additional files:\n')
      for f in addition:
        if f not in files:
          file.write('%s/%s\n' % (patch, f))
          files.append(f)

WritePatchFile('patch_files', ['drum'])
