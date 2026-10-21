# Contributing

GeneInIt is under active migration. Changes to biological selection or
post-annotation rules must include:

1. a versioned policy change,
2. an updated migration-matrix row,
3. a focused unit or adversarial-locus test, and
4. M2/M27 golden regression before a production tag.

Do not commit GeneMark keys, private containers, private genomes/RNA reads, or
absolute project paths.

Run unit tests with `python -m pytest tests/unit`. Synthetic fixtures are included
in the repository. The optional `tests/data/jinju_*.json` files are private local
regression fixtures and are excluded by `.gitignore`. Tests requiring unavailable
private fixtures are reported as skipped; synthetic tests still run. When those
files are present, the corresponding regression tests run normally. A public test
run therefore does not validate the private biological cases.
