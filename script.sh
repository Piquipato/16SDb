#!/bin/bash
CPUS=64
MEM=512
TIME_FIT="1-00:00:00"
PARTITION="gpu_h100_64C_128T_4TB_co_pi"
CONDA_ENV="propred"

for i in "$@"; do
    case $i in
        -c=*|--cpus=*)
        export CPUS="${i#*=}"
        shift # past argument=value
        ;;
        -m=*|--mem=*)
        export MEM="${i#*=}G"
        shift # past argument=value
        ;;
        -e=*|--conda-env=*)
        export CONDA_ENV="${i#*=}"
        shift # past argument=value
        ;;
        -t=*|--time=*)
        export TIME_FIT="${i#*=}"
        shift # past argument=value
        ;;
        -p=*|--partition=*)
        export PARTITION="${i#*=}"
        shift # past argument=value
        ;;
        *)
        # unknown option
        ;;
    esac
done

export SUBMIT_DIR="$(pwd)"
if [[ -e "$SUBMIT_DIR/work" && -d "$SUBMIT_DIR/work" ]]; then
    rm -rf "$SUBMIT_DIR/work"
fi
if [[ ! -e "$SUBMIT_DIR/work" ]]; then
    mkdir -p "$SUBMIT_DIR/work"
fi

export OUTPUT="$SUBMIT_DIR/work/database.fasta"
export LOGS="$SUBMIT_DIR/slurm-logs/fastas-%j.out"
sbatch << EOT
#!/bin/bash
#SBATCH --partition=$PARTITION
#SBATCH --job-name=fastas
#SBATCH --nodes=1
#SBATCH --cpus-per-task=$CPUS
#SBATCH --mem=$MEM
#SBATCH --time=$TIME_FIT
#SBATCH --output=$LOGS

# Script
ml purge && ml GCC/14.2.0 && ml GCCcore/14.2.0
cd $SUBMIT_DIR
conda run -n propred python createDatabase.py --output $OUTPUT
EOT