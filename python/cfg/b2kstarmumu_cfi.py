"""
Configuration File Initializer for B2KstarMuMu analysis

"""

import FWCore.ParameterSet.Config as cms


####### 2011, 5E32 1E33,1.4E33 ########
HLTBitNames_DoubleMu_v1 = cms.vstring(
	"HLT_DoubleMu3_Jpsi_v1", # 5E32, v4.2-v5.3
	"HLT_DoubleMu3_Jpsi_v2", # 5E32, v6.1-v6.2
	"HLT_Dimuon6p5_Jpsi_v1", #5E32, v8.1-v8.3
	"HLT_Dimuon6p5_Jpsi_Displaced_v1", #5E32, v8.1-v8.3
	"HLT_Dimuon6p5_Barrel_Jpsi_v1", #5E32, v8.1-v8.3
	"HLT_DoubleMu3_Quarkonium_v1", #5E32, v4.2-v5.3
	"HLT_DoubleMu3_Quarkonium_v2", #5E32, v6.1-v6.2
	"HLT_DoubleMu3_LowMass_v1", #5E32, v6.1-v6.2
	"HLT_Dimuon6p5_LowMass_v1", #5E32, v8.1-v8.3
	"HLT_Dimuon6p5_LowMass_Displaced_v1", #5E32, v8.1-v8.3
	"HLT_DoubleMu3_Bs_v1", #5E32, v4.2-v5.3
	"HLT_DoubleMu2_Bs_v1", #5E32, v6.1-v6.2
	"HLT_DoubleMu2_Bs_v2", #5E32, v8.1-v8.3
	"HLT_DoubleMu2_Bs_v3", #1E33, v1.2-v2.4
	"HLT_DoubleMu2_Bs_v4", #1E33, v2.5
	"HLT_DoubleMu2_Bs_v5", #1.4E33, v1.2
	"HLT_Dimuon0_Jpsi_v1", #1E33, v1.2-v2.4
	"HLT_Dimuon0_Jpsi_v2", #1E33, v2.5
	"HLT_Dimuon0_Jpsi_v3", #1.4E33, v1.2
	"HLT_Dimuon0_Jpsi_v5", #2E33, v1.1-v1.2
	"HLT_Dimuon0_Jpsi_v6", #3E33,  v1.1-v5.0
	"HLT_Dimuon0_Jpsi_v9", #5E33,  v1.4
	"HLT_Dimuon0_Jpsi_NoVertexing_v2", #2E33, v1.1-v1.2
	"HLT_Dimuon0_Jpsi_NoVertexing_v3", #3E33, v1.1-v5.0
	"HLT_Dimuon0_Jpsi_NoVertexing_v6", #5E33, v1.4
	"HLT_Dimuon0_Upsilon_v1", #1E33, v1.2-v2.4
	"HLT_Dimuon0_Upsilon_v2",#1E33, v2.5
	"HLT_Dimuon0_Upsilon_v3",#1.4E33, v1.2
	"HLT_Dimuon0_Upsilon_v5",#2E33, v1.1-v1.2
	"HLT_Dimuon0_Upsilon_v6",#3E33, v1.1-v5.0
	"HLT_Dimuon0_Upsilon_v9",#5E33, v1.4
	"HLT_Dimuon4_Bs_Barrel_v2",#1E33, v1.2-v1.3
	"HLT_Dimuon4_Bs_Barrel_v3", #1E33, v2.2-v2.4
	"HLT_Dimuon4_Bs_Barrel_v4", #1E33, v2.5
	"HLT_Dimuon4_Bs_Barrel_v5",#1.4E33, v1.2
	"HLT_Dimuon4_Bs_Barrel_v7",#2E33, v1.1-v1.2
	"HLT_DoubleMu4_Dimuon4_Bs_Barrel_v1", #3E33, v1.1-v5.0
	"HLT_DoubleMu4_Dimuon4_Bs_Barrel_v4",#5E33, v1.4
	"HLT_Dimuon5_Upsilon_Barrel_v1", #1E33, v1.2-v2.4
	"HLT_Dimuon5_Upsilon_Barrel_v2",#1E33, v2.5
	"HLT_Dimuon5_Upsilon_Barrel_v3",#1.4E33, v1.2
	"HLT_Dimuon5_Upsilon_Barrel_v5",#2E33, v1.1-v1.2
	"HLT_Dimuon7_Upsilon_Barrel_v1",#3E33, v1.1-v5.0
	"HLT_Dimuon7_Upsilon_Barrel_v4",#5E33, v1.4
	"HLT_Dimuon9_Upsilon_Barrel_v1",#3E33, v1.1-v5.0
	"HLT_Dimuon9_Upsilon_Barrel_v4",#5E33, v1.4
	"HLT_Dimuon6_Bs_v1",#1E33, v1.2-v1.3
	"HLT_Dimuon6_Bs_v2",#1E33, v2.2-v2.4
	"HLT_Dimuon6_Bs_v3",#1E33, v2.5
	"HLT_Dimuon6_Bs_v4",#1.4E33, v1.2
	"HLT_Dimuon6_Bs_v6",#2E33, v1.1-v1.2
	"HLT_DoubleMu4_Dimuon6_Bs_v1",#3E33, v1.1-v5.0         
	"HLT_DoubleMu4_Dimuon6_Bs_v4",#5E33, v1.4         
	"HLT_Dimuon7_LowMass_Displaced_v1",#1E33, v1.2-v1.3
	"HLT_Dimuon7_LowMass_Displaced_v2", #1E33, v2.2-v2.4
	"HLT_Dimuon7_LowMass_Displaced_v3",#1E33, v2.5
	"HLT_Dimuon7_LowMass_Displaced_v4",#1.4E33, v1.2
	"HLT_DoubleMu4_LowMass_Displaced_v2",#2E33, v1.1-v1.2
	"HLT_DoubleMu4p5_LowMass_Displaced_v1",#3E33, v1.1-v5.0
	"HLT_DoubleMu4p5_LowMass_Displaced_v4",#5E33, v1.4
	"HLT_DoubleMu5_LowMass_Displaced_v1",#3E33, v1.1-v5.0
	"HLT_DoubleMu5_LowMass_Displaced_v4",#5E33, v1.4
	"HLT_Dimuon7_Jpsi_Displaced_v1", #1E33, v1.2-v2.4
	"HLT_Dimuon7_Jpsi_Displaced_v2",#1E33, v2.5
	"HLT_Dimuon7_Jpsi_Displaced_v3",#1.4E33, v1.2
	"HLT_DoubleMu3p5_Jpsi_Displaced_v2", #2E33, v1.1-v1.2
	"HLT_DoubleMu4_Jpsi_Displaced_v1", #3E33, v1.1-v5.0
	"HLT_DoubleMu4_Jpsi_Displaced_v4", #5E33, v1.4
	"HLT_DoubleMu5_Jpsi_Displaced_v1", #3E33, v1.1-v3.1
	"HLT_DoubleMu5_Jpsi_Displaced_v2", #3E33, v4.0-v5.0
	"HLT_DoubleMu5_Jpsi_Displaced_v4", #5E33, v1.4
	"HLT_Dimuon7_Jpsi_X_Barrel_v1", #1E33, v1.2-v2.4
	"HLT_Dimuon7_Jpsi_X_Barrel_v2",#1E33, v2.5
	"HLT_Dimuon7_Jpsi_X_Barrel_v3",#1.4E33, v1.2
	"HLT_Dimuon7_Jpsi_X_Barrel_v5",#2E33, v1.1-v1.2
	"HLT_Dimuon7_PsiPrime_v1", #1E33, v1.2-v2.4
	"HLT_Dimuon7_PsiPrime_v2",#1E33, v2.5
	"HLT_Dimuon7_PsiPrime_v3",#1.4E33, v1.2
	"HLT_Dimuon7_PsiPrime_v5",#2E33, v1.1-v1.2
	"HLT_Dimuon9_PsiPrime_v1",#3E33, v1.1-v1.2
	"HLT_Dimuon9_PsiPrime_v4",#5E33, v1.4
	"HLT_Dimuon11_PsiPrime_v1",#3E33, v1.1-v1.2
	"HLT_Dimuon11_PsiPrime_v4",#5E33, v1.4
	"HLT_Dimuon10_Jpsi_Barrel_v1", #1E33, v1.2-v2.4
	"HLT_Dimuon10_Jpsi_Barrel_v2",#1E33, v2.5
	"HLT_Dimuon10_Jpsi_Barrel_v3",#1.4E33, v1.2
	"HLT_Dimuon10_Jpsi_Barrel_v5",#2E33, v1.1-v1.2
	"HLT_Dimuon10_Jpsi_Barrel_v6",#3E33, v1.1-v5.0
	"HLT_Dimuon10_Jpsi_Barrel_v9",#5E33, v1.4
	"HLT_Dimuon13_Jpsi_Barrel_v1",#3E33, v1.1-v5.0
	"HLT_Dimuon13_Jpsi_Barrel_v4",#5E33, v1.4
	"HLT_Dimuon0_Jpsi_Muon_v1", #1E33, v1.2-v1.3
	"HLT_Dimuon0_Jpsi_Muon_v2", #1E33, v2.2-v2.4
	"HLT_Dimuon0_Jpsi_Muon_v3",#1E33, v2.5
	"HLT_Dimuon0_Jpsi_Muon_v4",#1.4E33, v1.2
	"HLT_Dimuon0_Jpsi_Muon_v6", #2E33, v1.1-v1.2
	"HLT_Dimuon0_Jpsi_Muon_v7", #3E33, v1.1-v5.0
	"HLT_Dimuon0_Jpsi_Muon_v10", #5E33, v1.4
	"HLT_Dimuon0_Upsilon_Muon_v1",#1E33, v1.2-v1.3
	"HLT_Dimuon0_Upsilon_Muon_v2", #1E33, v2.2-v2.4
	"HLT_Dimuon0_Upsilon_Muon_v3",#1E33, v2.5
	"HLT_Dimuon0_Upsilon_Muon_v4",#1.4E33, v1.2
	"HLT_Dimuon0_Upsilon_Muon_v6",#2E33, v1.1-v1.2
	"HLT_Dimuon0_Upsilon_Muon_v7",#3E33, v1.1-v5.0
	"HLT_Dimuon0_Upsilon_Muon_v10",#5E33, v1.4
        )


# ONE FILTER NAME PER PATH    
HLTLastFilterNames_DoubleMu_v1 = cms.vstring(
	"hltDoubleMu3JpsiL3Filtered", #HLT_DoubleMu3_Jpsi_v1
	"hltDoubleMu3JpsiL3Filtered", #HLT_DoubleMu3_Jpsi_v2
	"hltDimuon6p5JpsiL3Filtered", #HLT_Dimuon6p5_Jpsi_v1
	"hltDimuon6p5JpsiDisplacedL3Filtered", #HLT_Dimuon6p5_Jpsi_Displaced_v1
	"hltDimuon6p5BarrelJpsiL3Filtered", #HLT_Dimuon6p5_Barrel_Jpsi_v1
	"hltDoubleMu3QuarkoniumL3Filtered", #HLT_DoubleMu3_Quarkonium_v1	      
	"hltDoubleMu3QuarkoniumL3Filtered", #HLT_DoubleMu3_Quarkonium_v2
	"hltDoubleMu3LowMassL3Filtered", #HLT_DoubleMu3_LowMass_v1
	"hltDimuon6p5LowMassL3Filtered", #HLT_Dimuon6p5_LowMass_v1
	"hltDimuon6p5LowMassL3FilteredDisplaced", #HLT_Dimuon6p5_LowMass_Displaced_v1
	"hltDoubleMu3BsL3Filtered", #HLT_DoubleMu3_Bs_v1
	"hltDoubleMu2BsL3Filtered", #HLT_DoubleMu2_Bs_v1
	"hltDoubleMu2BsL3Filtered", #HLT_DoubleMu2_Bs_v2
	"hltDoubleMu2BsL3Filtered", #HLT_DoubleMu2_Bs_v3
	"hltDoubleMu2BsL3Filtered", #HLT_DoubleMu2_Bs_v4
	"hltDoubleMu2BsL3Filtered", #HLT_DoubleMu2_Bs_v5
	"hltVertexmumuFilterJpsi", #HLT_Dimuon0_Jpsi_v1 was "hltJpsiL3Filtered"
	"hltVertexmumuFilterJpsi", #HLT_Dimuon0_Jpsi_v2 was "hltJpsiL3Filtered"
	"hltVertexmumuFilterJpsi", #HLT_Dimuon0_Jpsi_v3 was "hltJpsiL3Filtered"
	"hltVertexmumuFilterJpsi", #HLT_Dimuon0_Jpsi_v5
	"hltVertexmumuFilterJpsi", #HLT_Dimuon0_Jpsi_v6
	"hltVertexmumuFilterJpsi", #HLT_Dimuon0_Jpsi_v9
	"hltJpsiNoVertexingL3Filtered", #HLT_Dimuon0_Jpsi_NoVertexing_v2
	"hltJpsiNoVertexingL3Filtered", #HLT_Dimuon0_Jpsi_NoVertexing_v3
	"hltJpsiNoVertexingL3Filtered", #HLT_Dimuon0_Jpsi_NoVertexing_v6
	"hltVertexmumuFilterUpsilon", #HLT_Dimuon0_Upsilon_v1 was "hltUpsilonL3Filtered"
	"hltVertexmumuFilterUpsilon", #HLT_Dimuon0_Upsilon_v2 was "hltUpsilonL3Filtered"
	"hltVertexmumuFilterUpsilon", #HLT_Dimuon0_Upsilon_v3 was "hltUpsilonL3Filtered"
	"hltVertexmumuFilterUpsilon", #HLT_Dimuon0_Upsilon_v5
	"hltVertexmumuFilterUpsilon", #HLT_Dimuon0_Upsilon_v6
	"hltVertexmumuFilterUpsilon", #HLT_Dimuon0_Upsilon_v9
	"hltDoubleMu2BarrelBsL3Filtered", #HLT_Dimuon4_Bs_Barrel_v2
	"hltDoubleMu2BarrelBsL3Filtered", #HLT_Dimuon4_Bs_Barrel_v3
	"hltDoubleMu2BarrelBsL3Filtered", #HLT_Dimuon4_Bs_Barrel_v4
	"hltDoubleMu2BarrelBsL3Filtered", #HLT_Dimuon4_Bs_Barrel_v5
	"hltDoubleMu2BarrelBsL3Filtered", #HLT_Dimuon4_Bs_Barrel_v7
	"hltVertexmumuFilterBs4", #HLT_DoubleMu4_Dimuon4_Bs_Barrel_v1
	"hltVertexmumuFilterBs4", #HLT_DoubleMu4_Dimuon4_Bs_Barrel_v4
	"hltVertexmumuFilterUpsilonBarrel", #HLT_Dimuon5_Upsilon_Barrel_v1 was "hltBarrelUpsilonL3Filtered"
	"hltVertexmumuFilterUpsilonBarrel", #HLT_Dimuon5_Upsilon_Barrel_v2 was "hltBarrelUpsilonL3Filtered"
	"hltVertexmumuFilterUpsilonBarrel", #HLT_Dimuon5_Upsilon_Barrel_v3 was "hltBarrelUpsilonL3Filtered"
	"hltVertexmumuFilterUpsilonBarrel", #HLT_Dimuon5_Upsilon_Barrel_v5
	"hltVertexmumuFilterDimuon7UpsilonBarrel", #HLT_Dimuon7_Upsilon_Barrel_v1
	"hltVertexmumuFilterDimuon7UpsilonBarrel", #HLT_Dimuon7_Upsilon_Barrel_v4
	"hltVertexmumuFilterDimuon9UpsilonBarrel", #HLT_Dimuon9_Upsilon_Barrel_v1
	"hltVertexmumuFilterDimuon9UpsilonBarrel", #HLT_Dimuon9_Upsilon_Barrel_v4
	"hltDoubleMu2Dimuon6BsL3Filtered", #HLT_Dimuon6_Bs_v1
	"hltDoubleMu2Dimuon6BsL3Filtered", #HLT_Dimuon6_Bs_v2
	"hltDoubleMu2Dimuon6BsL3Filtered", #HLT_Dimuon6_Bs_v3
	"hltDoubleMu2Dimuon6BsL3Filtered", #HLT_Dimuon6_Bs_v4
	"hltDoubleMu2Dimuon6BsL3Filtered", #HLT_Dimuon6_Bs_v6
	"hltVertexmumuFilterBs6", #HLT_DoubleMu4_Dimuon6_Bs_v1 
	"hltVertexmumuFilterBs6", #HLT_DoubleMu4_Dimuon6_Bs_v4 
	"hltDisplacedmumuFilterLowMass", #HLT_Dimuon7_LowMass_Displaced_v1 was "hltLowMassDisplacedL3Filtered"
	"hltDisplacedmumuFilterLowMass", #HLT_Dimuon7_LowMass_Displaced_v2 was "hltLowMassDisplacedL3Filtered"
	"hltDisplacedmumuFilterLowMass", #HLT_Dimuon7_LowMass_Displaced_v3 was "hltLowMassDisplacedL3Filtered"
	"hltDisplacedmumuFilterLowMass", #HLT_Dimuon7_LowMass_Displaced_v4 was "hltLowMassDisplacedL3Filtered"
	"hltDisplacedmumuFilterLowMass", #HLT_DoubleMu4_LowMass_Displaced_v2
	"hltDisplacedmumuFilterDoubleMu4p5LowMass",#HLT_DoubleMu4p5_LowMass_Displaced_v1
	"hltDisplacedmumuFilterDoubleMu4p5LowMass",#HLT_DoubleMu4p5_LowMass_Displaced_v4
	"hltDisplacedmumuFilterDoubleMu5LowMass",#HLT_DoubleMu5_LowMass_Displaced_v1
	"hltDisplacedmumuFilterDoubleMu5LowMass",#HLT_DoubleMu5_LowMass_Displaced_v4
	"hltDisplacedmumuFilterJpsi", #HLT_Dimuon7_Jpsi_Displaced_v1 was "hltJpsiDisplacedL3Filtered"
	"hltDisplacedmumuFilterJpsi", #HLT_Dimuon7_Jpsi_Displaced_v2 was "hltJpsiDisplacedL3Filtered"
	"hltDisplacedmumuFilterJpsi", #HLT_Dimuon7_Jpsi_Displaced_v3 was "hltJpsiDisplacedL3Filtered"
	"hltDisplacedmumuFilterJpsi", #HLT_DoubleMu3p5_Jpsi_Displaced_v2
	"hltDisplacedmumuFilterDoubleMu4Jpsi", #HLT_DoubleMu4_Jpsi_Displaced_v1
	"hltDisplacedmumuFilterDoubleMu4Jpsi", #HLT_DoubleMu4_Jpsi_Displaced_v4
	"hltDisplacedmumuFilterDoubleMu5Jpsi", #HLT_DoubleMu5_Jpsi_Displaced_v1
	"hltDisplacedmumuFilterDoubleMu5Jpsi", #HLT_DoubleMu5_Jpsi_Displaced_v2
	"hltDisplacedmumuFilterDoubleMu5Jpsi", #HLT_DoubleMu5_Jpsi_Displaced_v4          
	"hltVertexmumuFilterJpsiXBarrel", #HLT_Dimuon7_Jpsi_X_Barrel_v1 was "hltJpsiXBarrelL3Filtered"
	"hltVertexmumuFilterJpsiXBarrel", #HLT_Dimuon7_Jpsi_X_Barrel_v2 was "hltJpsiXBarrelL3Filtered"
	"hltVertexmumuFilterJpsiXBarrel", #HLT_Dimuon7_Jpsi_X_Barrel_v3 was "hltJpsiXBarrelL3Filtered"
	"hltVertexmumuFilterJpsiXBarrel", #HLT_Dimuon7_Jpsi_X_Barrel_v5
	"hltVertexmumuFilterPsiPrime", #HLT_Dimuon7_PsiPrime_v1 was "hltPsiPrimeL3Filtered"
	"hltVertexmumuFilterPsiPrime", #HLT_Dimuon7_PsiPrime_v2 was "hltPsiPrimeL3Filtered"
	"hltVertexmumuFilterPsiPrime", #HLT_Dimuon7_PsiPrime_v3 was "hltPsiPrimeL3Filtered"
	"hltVertexmumuFilterPsiPrime", #HLT_Dimuon7_PsiPrime_v5
	"hltVertexmumuFilterDimuon9PsiPrime", #HLT_Dimuon9_PsiPrime_v1
	"hltVertexmumuFilterDimuon9PsiPrime", #HLT_Dimuon9_PsiPrime_v4
	"hltVertexmumuFilterDimuon11PsiPrime", #HLT_Dimuon11_PsiPrime_v1
	"hltVertexmumuFilterDimuon11PsiPrime", #HLT_Dimuon11_PsiPrime_v4
	"hltVertexmumuFilterJpsiBarrel", #HLT_Dimuon10_Jpsi_Barrel_v1 was "hltBarrelJpsiL3Filtered"
	"hltVertexmumuFilterJpsiBarrel", #HLT_Dimuon10_Jpsi_Barrel_v2 was "hltBarrelJpsiL3Filtered"
	"hltVertexmumuFilterJpsiBarrel", #HLT_Dimuon10_Jpsi_Barrel_v3 was "hltBarrelJpsiL3Filtered"
	"hltVertexmumuFilterJpsiBarrel", #HLT_Dimuon10_Jpsi_Barrel_v5
	"hltVertexmumuFilterDimuon10JpsiBarrel", #HLT_Dimuon10_Jpsi_Barrel_v6
	"hltVertexmumuFilterDimuon10JpsiBarrel", #HLT_Dimuon10_Jpsi_Barrel_v9
	"hltVertexmumuFilterDimuon13JpsiBarrel", #HLT_Dimuon13_Jpsi_Barrel_v1
	"hltVertexmumuFilterDimuon13JpsiBarrel", #HLT_Dimuon13_Jpsi_Barrel_v4
	"hltVertexmumuFilterJpsiMuon", #HLT_Dimuon0_Jpsi_Muon_v1 was "hltJpsiMuonL3Filtered"
	"hltVertexmumuFilterJpsiMuon", #HLT_Dimuon0_Jpsi_Muon_v2 was "hltJpsiMuonL3Filtered"
	"hltVertexmumuFilterJpsiMuon", #HLT_Dimuon0_Jpsi_Muon_v3 was "hltJpsiMuonL3Filtered"
	"hltVertexmumuFilterJpsiMuon", #HLT_Dimuon0_Jpsi_Muon_v4 was "hltJpsiMuonL3Filtered"
	"hltVertexmumuFilterJpsiMuon", #HLT_Dimuon0_Jpsi_Muon_v6
	"hltVertexmumuFilterJpsiMuon", #HLT_Dimuon0_Jpsi_Muon_v7
	"hltVertexmumuFilterJpsiMuon", #HLT_Dimuon0_Jpsi_Muon_v10
	"hltVertexmumuFilterUpsilonMuon", #HLT_Dimuon0_Upsilon_Muon_v1 was "hltUpsilonMuonL3Filtered"
	"hltVertexmumuFilterUpsilonMuon", #HLT_Dimuon0_Upsilon_Muon_v2 was "hltUpsilonMuonL3Filtered"
	"hltVertexmumuFilterUpsilonMuon", #HLT_Dimuon0_Upsilon_Muon_v3 was "hltUpsilonMuonL3Filtered"
	"hltVertexmumuFilterUpsilonMuon", #HLT_Dimuon0_Upsilon_Muon_v4 was "hltUpsilonMuonL3Filtered"
	"hltVertexmumuFilterUpsilonMuon", #HLT_Dimuon0_Upsilon_Muon_v6
	"hltVertexmumuFilterUpsilonMuon", #HLT_Dimuon0_Upsilon_Muon_v7
	"hltVertexmumuFilterUpsilonMuon", #HLT_Dimuon0_Upsilon_Muon_v10
	)


MuonSelection_v1 = cms.string("(isGlobalMuon || isTrackerMuon || (innerTrack.isNonnull && genParticleRef(0).isNonnull)) && abs(innerTrack.dxy)<4 && abs(innerTrack.dz)<35 && muonID('TrackerMuonArbitrated')")


