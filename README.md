# ZprimeTollMET

	cmsrel CMSSW_10_2_18
	cd CMSSW_10_2_18/src
	cmsenv
	git cms-init
	git cms-merge-topic yongbinfeng:DeepMET
	git cms-merge-topic cms-egamma:EgammaPostRecoTools
	git clone -b 102X git@github.com:hjkwon260/ZprimeTollMET.git
	scram b -j 8