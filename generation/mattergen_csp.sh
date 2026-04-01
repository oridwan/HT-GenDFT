#!/bin/bash
#SBATCH --job-name=mattergen_csp_train
#SBATCH --partition=GPU
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64G
#SBATCH --gres=gpu:2
#SBATCH --time=72:00:00
#SBATCH --output=mattergen_csp_%j.out
#SBATCH --error=mattergen_csp_%j.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=EMAIL_PLACEHOLDER

module purge
module load cuda/11.8
ulimit -s unlimited

source ~/.bashrc
conda activate mattersim

cd $HOME/SOFT/mattergen_test/

# W&B setup (offline mode)
export WANDB_PROJECT="crystal-generation"
export WANDB_JOB_TYPE="train"
export WANDB_NAME="csp_from_scratch"
export WANDB_MODE=offline

export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

echo "=========================================="
echo "MatterGen CSP Training (From Scratch)"
echo "=========================================="
echo "Config: CSP mode"
echo "Dataset: mp_20"
echo "Properties: None (pure CSP)"
echo "Strategy: DDP (2 GPUs)"
echo "Max epochs: 2000"
echo "=========================================="

# Train CSP model from scratch
# Use --config-name=csp to load the CSP configuration
mattergen-train \
    --config-name=csp \
    data_module=mp_20 \
    trainer.strategy=ddp_find_unused_parameters_true \
    trainer.devices=2 \
    trainer.num_nodes=1 \
    trainer.accelerator=gpu \
    trainer.precision=32 \
    trainer.accumulate_grad_batches=4 \
    trainer.max_epochs=2000 \
    data_module.num_workers.train=2 \
    data_module.num_workers.val=2

echo ""
echo "=========================================="
echo "CSP Training completed!"
echo "=========================================="

# Find the latest checkpoint
if [ -d "outputs/singlerun" ]; then
    echo ""
    echo "Checkpoint location:"
    LATEST_RUN=$(find outputs/singlerun -type d -name "20*" -maxdepth 2 | sort | tail -1)
    if [ -n "$LATEST_RUN" ]; then
        echo "  $LATEST_RUN"
        
        # Find W&B run directory
        WANDB_DIR=$(find "$LATEST_RUN" -type d -name "offline-run-*" 2>/dev/null | head -1)
        if [ -n "$WANDB_DIR" ]; then
            PARENT_DIR=$(dirname "$WANDB_DIR")
            CKPT_DIR="$PARENT_DIR/../crystal-generation"
            
            if [ -d "$CKPT_DIR" ]; then
                CKPT_SUBDIR=$(find "$CKPT_DIR" -type d -name "checkpoints" | head -1)
                if [ -n "$CKPT_SUBDIR" ]; then
                    echo ""
                    echo "Checkpoints:"
                    ls -lh "$CKPT_SUBDIR"
                    echo ""
                    echo "Best checkpoint path for generate_ternary_csp.sh:"
                    echo "  CSP_MODEL=\"$LATEST_RUN\""
                fi
            fi
        fi
    fi
fi

echo ""
echo "=========================================="
echo "Next steps:"
echo "  1. Update CSP_MODEL path in ternary_electride/generate_ternary_csp.sh"
echo "  2. Run: sbatch ternary_electride/generate_ternary_csp.sh"
echo "=========================================="

