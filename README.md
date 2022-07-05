# ZprimeTollMET

	cmsrel CMSSW_10_6_19_patch3
	cd CMSSW_10_6_19_patch3/src
	cmsenv
	git cms-init
	# temporary
	# git cms-merge-topic yongbinfeng:DeepMET
	# git cms-merge-topic cms-egamma:EgammaPostRecoTools
	git clone -b 102X git@github.com:hjkwon260/ZprimeTollMET.git
	scram b -j 8
