"""
Module for CRAB operations

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"

import sys
import atr
from tls import * 


def main(args):
    if args[0] == 'cfg':
        cfg(args[1:])
    elif args[0] == 'gen':
        gen(args[1:])
    else:
        raise NameError(args)

def cfg(args):
    datatype = args[0]
    label = args[1]

    if datatype == 'dt' and label in ['5ifbv2.3', '5ifbv2.4', '5ifbv2.5',
                                      '5ifbv2.6']:
        cfg = get_crab_cfg_template(label)
        datasets = atr.datasets(datatype, label)
        for dataset in datasets:
            filename = 'crab_%s.cfg' % dataset 
            key_val_dict = datasets[dataset]
            for k, v in key_val_dict.items():
                cfg.set_key_value(k, v)

            user_remote_dir = '/user/x/xshi/afb/dat/srm/%s/%s/%s' %(
                datatype, label, dataset)
            cfg.set_key_value('user_remote_dir', user_remote_dir)
            castor_dir = '/castor/cern.ch' + user_remote_dir
            check_and_join(castor_dir, '')
            cfg.output(filename, verbose=1)

    else:
        sys.stdout.write('name ui_working_dir in the cfg file! \n')
        raise NameError(args)


def gen(args):
    sample = args[0]
    label = args[1]

    cfg = get_crab_cfg_template(label)

    print cfg
    sys.exit()
    import ConfigParser

    config = ConfigParser.RawConfigParser()

    config.add_section('Section1')
    config.set('Section1', 'int', '15')
    config.set('Section1', 'bool', 'true')
    config.set('Section1', 'float', '3.1415')
    config.set('Section1', 'baz', 'fun')
    config.set('Section1', 'bar', 'Python')
    config.set('Section1', 'foo', '%(bar)s is %(baz)s!')

    dir(config)
    sys.exit()
    # Writing our configuration file to 'example.cfg'
    with open('example.cfg', 'wb') as configfile:
        config.write(configfile)

    sys.exit()    
        #print cfg
    cfg.output()

    if label == 'GEN_1M_v1.0':
        pass 
    sys.exit()
    
