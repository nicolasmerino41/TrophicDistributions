# Thermal master-data dictionary

## `thermal_traits_long.csv`

- `species_alias`, `accepted_species`: harmonized binomial names before and
  after GloBI-supported alias alignment.
- `metric`, `trait_value`: thermal metric and value in degrees Celsius.
- `trait_source`, `source_record_id`: provenance of the measurement.
- taxonomy and habitat fields: information retained from the source database.

## `thermal_species_master.csv`

- one row per species appearing in a trait source or retained interaction;
- `trait_value_*`, `trait_sd_*`, `trait_n_*`, and
  `selected_trait_source_*`: analysis values and provenance by metric;
- `source_trait_mean_*`, `source_trait_sd_*`, and `source_trait_n_*`:
  source-specific summaries retained for auditing;
- `recorded_consumer_degree`: unique non-self resources documented by GloBI or
  TetraEU;
- `recorded_resource_consumers`: number of documented consumers;
- `taxon_status`: animal, non-animal, or unresolved from the available
  taxonomic evidence;
- `globi_taxon_paths`, `globi_taxon_ids`, and `tetraeu_taxon_ids`: interaction
  database identifiers and taxonomic evidence.

## `thermal_recorded_interactions.csv`

- `consumer`, `resource`: harmonized species pair;
- `GloBI`, `TetraEU`, `interaction_sources`: database provenance;
- `interaction_types`, study and reported-name fields: retained evidence;
- `consumer_taxon_status`, `consumer_taxon_evidence`: validation of the
  consumer role. Only verified animals are retained as consumers; any taxon may
  remain as a resource;
- `self_link`: whether consumer and resource resolve to the same species.

## `thermal_interactions_master.csv`

The recorded-interaction fields plus both partners' selected thermal values,
sources, measurement counts, and metric-specific eligibility flags. Self-links
are retained in the recorded table but excluded from every thermal analysis.

## `thermal_consumer_coverage.csv`

- `recorded_prey_species`: unique documented non-self resources;
- `prey_with_trait`: recorded resources with the same thermal metric as the
  consumer;
- `prey_trait_coverage`: `prey_with_trait / recorded_prey_species`.

## `thermal_consumer_taxon_audit.csv`

- one row per thermal-trait species considered for a GloBI consumer query;
- `consumer_taxon_status`: animal, non-animal, or unresolved;
- `queried_as_consumer` and `exclusion_reason`: the reproducible decision and
  its reason. Only verified animals are queried as consumers.
