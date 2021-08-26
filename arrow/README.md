# Table conversion

```bash
analysis-runner --dataset seqr --access-level standard --output-dir seqr-conversion_$(date +"%Y-%m-%d_%H-%M-%S") --description "seqr table conversion" main.py --input=gs://path/to/annotated_input.mt
```
