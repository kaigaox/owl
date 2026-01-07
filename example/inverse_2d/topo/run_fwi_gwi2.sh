
p=param_dtw2.rb
cp -rp param_fwi.rb $p
echo 'misfit_type = dtw-vector ' >> $p
echo 'dir_working = test_dtw_vector ' >> $p
echo 'adj_nt = 1000 ' >> $p
echo 'dtw_form = syn-obs' >> $p
echo 'tlag_max = 0.25' >> $p
echo 'dtw_smooth_median = 1' >> $p
echo 'dtw_smooth_gaussian = 3' >> $p
echo 'dtw_rinst = 1.0 ' >> $p
echo 'dtw_rcuml = 1.0 ' >> $p
echo 'dtw_loss = l0.5 ' >> $p
echo 'jumpout_factor = 1.05 ' >> $p
$HOME/intel/mpi/bin/mpirun -np 60 $HOME/bin/owl_fwi2 $p
