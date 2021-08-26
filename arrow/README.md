# Loading pipeline

Usually the seqr loading pipeline annotates a Hail MatrixTable and subsequently loads that into Elasticsearch.

To instead convert the MatrixTable to the Arrow format, run the following table conversion in Hail Batch:

```bash
analysis-runner --dataset seqr --access-level standard --output-dir seqr_table_conversion/$(date +"%Y-%m-%d_%H-%M-%S") --description "seqr table conversion" main.py --input=gs://path/to/annotated_input.mt
```
