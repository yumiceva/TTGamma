#./ttgamma data muon inputskim outdir=./histos/ >& histos/histos_data.log &

./ttgamma data_mu_a muon inputskim outdir=./histos/ >& histos/histos_data_mu_a.log &
./ttgamma data_mu_b muon inputskim outdir=./histos/ >& histos/histos_data_mu_b.log &
./ttgamma data_mu_c muon inputskim outdir=./histos/ >& histos/histos_data_mu_c.log &
./ttgamma data_mu_d muon inputskim outdir=./histos/ >& histos/histos_data_mu_d.log &

./ttgamma t_s muon inputskim outdir=./histos/ >& histos/histos_t_s.log &
./ttgamma t_t muon inputskim outdir=./histos/ >& histos/histos_t_t.log &
./ttgamma t_tW muon inputskim outdir=./histos/ >& histos/histos_t_tW.log &
sleep 1m
./ttgamma tbar_s muon inputskim outdir=./histos/ >& histos/histos_tbar_s.log &
./ttgamma tbar_t muon inputskim outdir=./histos/ >& histos/histos_tbar_t.log &
./ttgamma tbar_tW muon inputskim outdir=./histos/ >& histos/histos_tbar_tW.log &
sleep 1m
./ttgamma DYJetsToLL muon inputskim outdir=./histos/ >& histos/histos_DYJetsToLL.log &
./ttgamma Wjets muon inputskim outdir=./histos/ >& histos/histos_Wjets.log &
sleep 1m
./ttgamma ttg muon inputskim outdir=./histos/ >& histos/histos_ttg.log &
./ttgamma ttgWhizard muon inputskim outdir=./histos/ >& histos/histos_ttgWhizard.log &
./ttgamma ttW muon inputskim outdir=./histos/ >& histos/histos_ttW.log &
./ttgamma ttZ muon inputskim outdir=./histos/ >& histos/histos_ttZ.log &
sleep 1m
./ttgamma ttjets_1l muon inputskim outdir=./histos/ >& histos/histos_ttjets_1l.log &
./ttgamma ttjets_2l muon inputskim outdir=./histos/ >& histos/histos_ttjets_2l.log &
./ttgamma ttjets_0l muon inputskim outdir=./histos/ >& histos/histos_ttjets_0l.log &

./ttgamma ttjets_1l muon inputskim onlyphotons outdir=./histos/ >& histos/histos_ttjets_1l_g.log &
./ttgamma ttjets_2l muon inputskim onlyphotons outdir=./histos/ >& histos/histos_ttjets_2l_g.log &
./ttgamma ttjets_0l muon inputskim onlyphotons outdir=./histos/ >& histos/histos_ttjets_0l_g.log &