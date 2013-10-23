"""
Attributes for the SRM files.

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"


def srmpath(label, datatype='dt'):
    if label in ['5ifbv2.3']:
        pbs = 'srm://srm.ihep.ac.cn:8443/srm/managerv2?SFN='
        fbs = '/afs/cern.ch/user/z/zhlinl/public/forXin/file_26Dec/'
        fs = [fbs + 'promptReco_V4' ]

        p = [ pbs+ '/pnfs/ihep.ac.cn/data/cms/store/user/zhlinl/B02Kstarmumu/PAT/Tree_26Dec/May10ReReco-v1/zhlinl/MuOnia/Run2011A-May10ReReco-v1-TTree-b2MuMuX-V5/f991a63cdabf1a4b8902294cffd8d647/',
              pbs+  '/pnfs/ihep.ac.cn/data/cms/store/user/zhlinl/B02Kstarmumu/PAT/Tree_26Dec/Run2011B-PromptReco-v1/zhlinl/MuOnia/Run2011B-PromptReco-v1-TTree-b2MuMuX-V5/2fd4dde7c1ebb37c2c287e61c3c3cc2c/', 

              pbs + '/pnfs/ihep.ac.cn/data/cms/store/user/zhlinl/B02Kstarmumu/PAT/Tree_26Dec/Run2011A-PromptReco-v4/zhlinl/MuOnia/Run2011A-PromptReco-v4-TTree-b2MuMuX-V5/2fd4dde7c1ebb37c2c287e61c3c3cc2c/', 

              pbs + '/pnfs/ihep.ac.cn/data/cms/store/user/zhlinl/B02Kstarmumu/PAT/Tree_26Dec/Run2011A-PromptReco-v5/zhlinl/MuOnia/Run2011A-PromptReco-v5-TTree-b2MuMuX-V5/2fd4dde7c1ebb37c2c287e61c3c3cc2c/', 

              pbs + '/pnfs/ihep.ac.cn/data/cms/store/user/zhlinl/B02Kstarmumu/PAT/Tree_26Dec/Run2011A-PromptReco-v6/zhlinl/MuOnia/Run2011A-PromptReco-v6-TTree-b2MuMuX-V5/2fd4dde7c1ebb37c2c287e61c3c3cc2c/'
              ]

    elif label in ['B2KstarMuMu/GEN_1M_v1', 'B2KstarMuMu/HLT_1M_v1',
                   'B2KstarMuMu/HLT_1M_v2', 'B2KstarMuMu/RECO_1M_v2']:
        p = ['/castor/cern.ch/user/x/xshi/afb/dat/srm/%s/%s' % (
            datatype, label)]

    elif label in ['B2KstarMuMu/RECO_1M_v2.4',
                   'B2KstarMuMu/GEN_SIM_10M_v2']:
        p = ['/store/user/xshi/afb/dat/srm//%s/%s' % (
            datatype, label)]

    elif label in ['B2KstarMuMu/RECO_10M_v2.1']:
        p = ['/store/user/xshi/B2KstarMuMu_GEN_SIM_10M_v2/B2KstarMuMu_RECO_10M_v2/74adfdc7814766ad34dc5589b8cb6432/']
        
    elif label in ['B2KstarMuMu/RECO_100M_v1.1']:
        p = ['/store/user/xshi/B2KstarMuMu_GEN_SIM_HLT_100M_v1/B2KstarMuMu_RECO_100M_v1/74adfdc7814766ad34dc5589b8cb6432/']
        
    else:
        raise NameError(label)
    
    return p

