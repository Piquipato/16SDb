#!/bin/bash
CPUS=32
MEM=512
TIME_FIT="1-00:00:00"
PARTITION="hmem_128C_256T_2TB"
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
        echo "Unknown option"
        exit 1
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
export CMDL="$(dirname $LOGS)/fastas-logs-$(date +"%Y_%m_%d_%I_%M_%S_%p").txt"
sbatch << EOT
#!/bin/bash
#SBATCH --partition=$PARTITION
#SBATCH --job-name=fastas
#SBATCH --nodes=1
#SBATCH --cpus-per-task=$CPUS
#SBATCH --mem="$MEM"G
#SBATCH --time=$TIME_FIT
#SBATCH --output=$LOGS

# Script
cd $SUBMIT_DIR
conda run -n propred python -u createDatabase.py --output $OUTPUT --n_jobs $CPUS >> $CMDL
EOT