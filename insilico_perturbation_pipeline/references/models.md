# Available State Models for In-Silico Perturbation

## Model Registry

| Model | Size | Config | Training Data | Status |
|-------|------|--------|---------------|--------|
| state_sm | Small | `state_sm` | Replogle (40k steps) | Tested |
| state_base | Base | `state` | Full dataset (10k steps) | Tested |

## Model Paths

### State SM (Small)
- **Dir**: `/nfs/turbo/umms-drjieliu/usr/zheyuz/state_with_esm/output_state_sm_baseline/state_sm_40k_20251202_001509`
- **Checkpoint**: `checkpoints/step=step=18000-val_loss=val_loss=1.7692.ckpt`
- **Notes**: Trained on Replogle data. Faster inference. Good for initial runs.

### State Base
- **Dir**: `/nfs/turbo/umms-drjieliu/usr/zheyuz/state_with_esm/output_state_base_full_10k/state_base_full_10k_20251206_010512`
- **Checkpoint**: `checkpoints/step=step=3000-val_loss=val_loss=1.5983.ckpt`
- **Notes**: Larger model, trained on full dataset. Better predictions but slower. Needs H100 GPU (use `--partition=drjieliu-h100`). Set `--time=20:00:00` for SLURM.

## Environment

- **Conda env**: `/nfs/turbo/umms-drjieliu/usr/zheyuz/miniforge-pypy3/envs/state_env`
- **Activate**: `source activate /nfs/turbo/umms-drjieliu/usr/zheyuz/miniforge-pypy3/envs/state_env`

## Gene Reference

- **18,080 model genes**: `/nfs/turbo/umms-drjieliu/usr/zheyuz/state_with_esm/competition_support_set/competition_support_set/gene_names.csv`
- Format: one gene symbol per line, no header

## Adding New Models

When a new model is trained, add an entry above with:
1. Model directory path
2. Best checkpoint path (lowest val_loss)
3. Training details (data, steps, architecture)
4. GPU requirements and estimated inference time
