./ttgamma data muon inputskim >& histos_data.log &
./ttgamma t_s muon inputskim >& histos_t_s.log &
./ttgamma t_t muon inputskim >& histos_t_t.log &
./ttgamma t_tW muon inputskim >& histos_t_tW.log &
sleep 1m
./ttgamma tbar_s muon inputskim >& histos_tbar_s.log &
./ttgamma tbar_t muon inputskim >& histos_tbar_t.log &
./ttgamma tbar_tW muon inputskim >& histos_tbar_tW.log &
sleep 1m
./ttgamma DYJetsToLL muon inputskim >& histos_DYJetsToLL.log &
./ttgamma Wjets muon inputskim >& histos_Wjets.log &
sleep 2m
./ttgamma diphoton_box_10to25 muon inputskim >& histos_diphoton_box_10to25.log &
./ttgamma diphoton_box_25to250 muon inputskim >& histos_diphoton_box_25to250.log &
./ttgamma diphoton_box_25toInf muon inputskim >& histos_diphoton_box_25toInf.log &
sleep 1m
./ttgamma ttjets_1l muon inputskim >& histos_ttjets_1l.log &
./ttgamma ttjets_2l muon inputskim >& histos_ttjets_2l.log &

