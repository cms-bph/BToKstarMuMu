# Package for B To K* Mu Mu Analysis

* NTuple production

  - MC KstarJPsi with BToKstarMuMu v2

    1. Get a small sample, edit the cfg/copy_cfg.py and run:

       : cd python 
       : cmsRun cfg/copy_cfg.py

    2. Generate the nutple command:
       : ./afb.py ntp mc BuToKstarJPsi-8TeV-2p5E7 -t 

       Here the "-t" means "test". The output is:
       : cmsRun btokstarmumu_MC.py
       : Please test with the above command.
       
       Make sure you can run the interactive cmsRun first.

    3. Create the crab cfg:
       : ./afb.py ntp mc BuToKstarJPsi-8TeV-2p5E7 -grid     
       
       Output:
       : Please create grid jobs with command: 
       : crab -cfg crab_BuToKstarJPsi_8TeV_2p5E7.cfg -create

       Then you can submit by: 
       : crab -c crab_BuToKstarJPsi_8TeV_2p5E7 -submit all 
       
       Note: make the job numbers less than 500 for each submission.
       Break down if necessary: 1-500, 501-1000, etc. 



 




 
