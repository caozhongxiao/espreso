#!/bin/bash

#module load cuda


#module load impi/4.1.1.036
#module load openmpi/1.8.1-icc
module load bullxmpi/bullxmpi_1.2.4.1


#module load intel/14.0.1
module load intel/15.2.164

#. /apps/intel/BETA2015/haswell.sh

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./libs

export OMP_NUM_THREADS=16

export MKL_PARDISO_OOC_MAX_CORE_SIZE=3500
export MKL_PARDISO_OOC_MAX_SWAP_SIZE=2000

log_file=output_8_nodes.log

dom_size=20      #20
clusters=3       #2
domains_p_c=5    #7

#         HEXA8 TETRA4 TETRA10
el_type=(   0     1      2    )

#               OM OK OK 
#               0   1   2   3   4   5   6   7   8   9
dom_size=(      10   8   3  10   6  14  16  15  10  10 )

clustt_size_x=( 2   8   2   8   16   9   3   3   3   5 )
clustt_size_y=( 2   8   2   8   16   9   3   4   3   5 )
clustt_size_z=( 2   8   2   8   16   9   8   8   3   5 ) 

clusters_x=(    1   1   2   1   1   2   8   8   1   2 )
clusters_y=(    1   1   2   1   1   2   8   8   1   2 )
clusters_z=(    1   1   2   1   1   2   4   4   1   2 )

corners=(       0   0   0   0   3   2   3   3   3   3 )

for i in 0 # 0 1 2 3
do
   d=${dom_size[${i}]}
   c=${corners[${i}]}   

   x=${clustt_size_x[${i}]}
   y=${clustt_size_y[${i}]}
   z=${clustt_size_z[${i}]}
  
   X=${clusters_x[${i}]}
   Y=${clusters_y[${i}]}
   Z=${clusters_z[${i}]}


   log_file=LOG-$X:$Y:$Z-$x:$y:$z-$d:$d:$d-$c:$c:$c.log

   date | tee $log_file
 
   echo "Config: dom_size=  $d | cluster_size = $x:$y:$z | clusters = $X:$Y:$Z  " 

   date | tee -a $log_file

   # ./espreso ${el_type[0]} ${x} ${y} ${z} ${d} ${d} ${d} | tee -a $log_file
   mpirun -bind-to-none -n $(( X * Y * Z ))  ./espreso ${el_type[0]} ${x} ${y} ${z} ${d} ${d} ${d} | tee -a $log_file
done

