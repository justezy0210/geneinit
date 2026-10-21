"""Reuse trained AUGUSTUS to add evidence-guided candidates in masked regions."""
from pathlib import Path
import json
import sys

sys.path.insert(0, str(snakemake.params.repo_root))
from geneinit.contracts.deployment import apptainer_command, load_predictor_deployment
from geneinit.evidence.components import load_junction_aggregate, load_protein_chains, load_repeats, load_stringtie_transcripts, load_trace
from geneinit.models.candidate import parse_predictor_gtf
from geneinit.predictor.common import sha256_file, write_completion
from geneinit.predictor.commands import augustus_prediction_command
from geneinit.predictor.masking import METHOD, discover_targets, masked_fasta, prepare_windows, collect_predictions
from geneinit.predictor.parallel import run_prediction_jobs

work=Path(snakemake.output.work);work.mkdir(parents=True,exist_ok=True)
Path(snakemake.log[0]).parent.mkdir(parents=True,exist_ok=True)
deployment=load_predictor_deployment(dict(snakemake.config))
genome=masked_fasta(Path(snakemake.input.genome))
audit_genome={key:seq.upper() for key,seq in genome.items()}
models=parse_predictor_gtf(Path(snakemake.input.augustus),'augustus',audit_genome)
models+=parse_predictor_gtf(Path(snakemake.input.genemark),'genemark_etp',audit_genome)
del audit_genome
transcripts=load_stringtie_transcripts(Path(snakemake.input.transcripts))
chains=load_protein_chains(Path(snakemake.input.direct_chains),'direct')
chains+=load_protein_chains(Path(snakemake.input.nonhc_chains),'nonhc',trace=load_trace(Path(snakemake.input.trace)))
targets=discover_targets(models,transcripts,chains,load_repeats(Path(snakemake.input.repeats)),load_junction_aggregate(Path(snakemake.input.junctions)))
windows,shards=prepare_windows(targets,genome,Path(snakemake.input.hints),work)
Path(snakemake.output.targets).write_text(json.dumps({'method':METHOD,'targets':windows},indent=2)+'\n')
def command(shard):
    inner=augustus_prediction_command(species=str(snakemake.params.species),genome=shard.fasta,hints=shard.hints,
        extrinsic_config=Path(snakemake.input.extrinsic),augustus_config=Path(snakemake.input.augustus_config),output=shard.output)
    return apptainer_command(deployment,inner,(work,Path(snakemake.input.extrinsic),Path(snakemake.input.augustus_config)),pwd=work)
maximum=run_prediction_jobs(shards,int(snakemake.threads),command,Path(snakemake.log[0])) if shards else 0
if not shards:Path(snakemake.log[0]).write_text('No evidence-supported masking conflicts found.\n')
counts=collect_predictions(windows,shards,genome,Path(snakemake.output.gtf))
Path(snakemake.output.summary).write_text(json.dumps({'status':'PASS','method':METHOD,'targets':len(windows),
    'maximum_parallel_augustus':maximum,'retrained':False,'independent_predictor_credit':False,**counts},indent=2)+'\n')
write_completion(Path(snakemake.output.complete),[('status','PASS'),('method',METHOD),('targets',len(windows)),
    ('candidate_models',counts.get('candidate_models',0)),('gtf_sha256',sha256_file(Path(snakemake.output.gtf))),
    ('hints_sha256',sha256_file(Path(snakemake.input.hints))),('extrinsic_sha256',sha256_file(Path(snakemake.input.extrinsic))),
    ('container_sha256',deployment.sif_sha256)])
