# assignment-7
module load craype-accel-nvidia80
srun -p gpu --gres=gpu:1 --ntasks=1 --time=00:05:00 --mem=40G --reservation=p_covidpre_121 ./laplace

