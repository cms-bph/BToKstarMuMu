"""
Module for Clean Data files 

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"

import sys 
from tls import * 
import atr 
import json


def main(args):
    if args[0] == 'eos':
        eos_rm(args[1])
    elif args[0] == 'data' or args[0] == 'mc':
        label = args[1]
        eosdir = atr.ntp.grid_path(label)
        quick_rm(eosdir)
    else:
        raise NameError(args)

def eos_rm(dir):
    dbname = dir.replace('/', '_')
    dbpath = os.path.join(atr.datpath, 'eos', 'db')
    dbfile = check_and_join(dbpath, dbname)

    if not os.access(dbfile, os.F_OK):
        sys.stdout.write('Create the db file ... \n')
        cmd = 'cmsLs %s' % dir
        output = proc_cmd(cmd)
        files = parse_ls(output)
        db = open(dbfile, 'w')    
        json.dump(files, db)
        db.close()
    else:
        db = open(dbfile)
        files = json.load(db)
        db.close()
 
    sys.stdout.write('Cleaning %s files... ' %len(files))
    sys.stdout.flush()
    
    n = 0 
    for f in files:
        n += 1 
        sys.stdout.write('%s  ' %n)
        sys.stdout.flush()
        cmd = 'cmsRm %s ' % f
        proc_cmd(cmd, procdir=dir)

    sys.stdout.write(' all done.\n')

    
def parse_ls(content, ignore=None, poz=None):
    names = []
    content = content.strip()
    lines = content.split('\n')
    for l in lines:
        if ignore != None and ignore in l:
                continue
        name = l.split()[-1]
        names.append(name)
    return names

def quick_rm(dir):
    cmd = 'rm -rf %s' % dir
 
    sys.stdout.write('Have you checked the output been merged correctly? \n')
    s = raw_input('Do you want to remove all of raw output files? (yes or no) ')
    if s != 'yes':
        return

    proc_cmd(cmd)
    sys.stdout.write(' done.\n')     
