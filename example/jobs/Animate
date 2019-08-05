#!/bin/bash
#SBATCH -t 01:00:00
#SBATCH --output=../stds/synth_%j.out
#SBATCH --error=../stds/synth_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=joeri.vanengelen@deltares.nl
#SBATCH -p short

set -x
echo $TMPDIR

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/visualization/llvm/lib
export PATH=$HOME/visualization/llvm/bin:$PATH

export PATH=$HOME/visualization/paraview/bin:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/visualization/paraview/lib

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/visualization/mesa/lib
export PATH=$HOME/visualization/mesa/bin:$PATH

export PYTHONPATH=$PYTHONPATH:$HOME/visualization/paraview/lib/paraview-5.4/site-packages

#outfol=$HOME/synthdelta/results/synth_SD_i202_m24_6577661_cont
outfol=$1

regex="([0-9]{7})"
if [[ ${outfol} =~ $regex ]]
then
    name="${BASH_REMATCH[1]}"
else
    name="results"
fi


mkdir -m 755 ${outfol}/screenshots

pvbatch --use-offscreen-rendering $HOME/synthdelta/delta_aquifer/post/paraview_macros/results_animate.py ${outfol}

#Create animation from snapshots
cd ${outfol}/screenshots

module load surf-devel
module load 2019
module load FFmpeg

filters="fps=5,scale=1000:-1:flags=lanczos"
ffmpeg -i results.%04d.png -vf "$filters,palettegen" -y palette.png
ffmpeg -framerate 5 -i results.%04d.png -i palette.png -lavfi "$filters [x]; [x][1:v] paletteuse" -r 30 -y "${name}.gif"

