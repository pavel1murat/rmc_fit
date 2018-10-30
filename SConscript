#!/usr/bin/env python

import os, re, string, subprocess

Import('env')
# Import('mu2e_helper')
#------------------------------------------------------------------------------
# last two components of the path. Ex: /not/this/but/THIS/AND_THIS
#                                      "AND_THIS" is usually "src"
#------------------------------------------------------------------------------
# x = subprocess.call('scripts/build_config',shell=True)
env['CXXFLAGS'].append('-I'+os.environ['MU2E_BASE_RELEASE']+'/include');
#------------------------------------------------------------------------------
# done
#------------------------------------------------------------------------------
from stntuple_helper import *

env.Append(BUILDERS = {'StntupleCodegen'  : stntuple_codegen})
env.Append(BUILDERS = {'StntupleRootCint' : stntuple_rootcint})

Export('stntuple_helper')
