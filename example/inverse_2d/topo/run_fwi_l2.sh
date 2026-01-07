
p=param_waveform.rb
cp -rp param_fwi.rb $p
echo 'misfit_type = waveform ' >> $p
echo 'dir_working = test_waveform ' >> $p
echo 'jumpout_factor = 1~40:1, 41:1.05' >> $p
$HOME/intel/mpi/bin/mpirun -np 60 $HOME/bin/owl_fwi2 $p
