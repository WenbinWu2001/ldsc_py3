# Resources and data availability

Last updated on: 2026-10-29

[Audit overview](README.md) · [Workflow counts and membership](workflow-index.md) · [Verification](verification.md)

## Main artifacts requested by users

The recurring request is for compatible reference resources for a particular ancestry, genome build and analysis goal.
Users ask both for individual missing files and for confirmation that files obtained from different releases or hosts belong together.
This summary draws on the issues reviewed in the 2026-10-28 audit; it does not establish that every requested combination is available or scientifically validated.

| Artifact family | What users are requesting | Example issues |
| --- | --- | --- |
| **PLINK reference panels** | Downloadable genotype panels, especially EUR Phase 3; ancestry/subpopulation matches, including CHB versus broader EAS; and the sample selection and SNP QC needed to reproduce a distributed panel. | [#414](#issue-414), [#479](#issue-479), [#495](#issue-495) |
| **Precomputed LD scores and regression weights** | Ordinary EUR/EAS scores for SNP heritability and genetic correlation, weights for non-European analyses, and resources beyond EUR/EAS; users need a mapping from old archive names to compatible replacements. | [#376](#issue-376), [#382](#issue-382), [#429](#issue-429), [#499](#issue-499), [#503](#issue-503) |
| **Complete hg19/GRCh37 and hg38/GRCh38 resource sets** | Matching LD scores, weights, baseline annotations, cell-type queries and HM3 coordinates/alleles for the chosen genome build; notably a complete EUR hg38 set, including baseline v1.2 and Cahoy resources. | [#264](#issue-264), [#331](#issue-331), [#357](#issue-357), [#470](#issue-470), [#505](#issue-505) |
| **Baseline annotation files and model definitions** | Versioned baseline/baseline-LD resources, original BED or full `.annot` files, corresponding LD scores, and a manifest explaining annotation names and biological meaning. | [#224](#issue-224), [#233](workflow.md#issue-233), [#479](#issue-479), [#505](#issue-505) |
| **Cell-type/tissue query annotations and LD scores** | Cahoy and other legacy cell-type collections, CNS examples, and GTEx brain versus multi-tissue annotations; users also need the matching baseline and control definitions. | [#106](#issue-106), [#201](#issue-201), [#248](#issue-248), [#412](#issue-412), [#450](#issue-450) |
| **HM3 restriction and allele-matching lists** | Files such as `w_hm3.snplist` and `hapmap3_snps.tgz`, hg38 coordinates/alleles, and the actual SNP lists used with EAS or AFR resources; these requests concern matching and provenance, not just a download basename. | [#235](#issue-235), [#322](#issue-322), [#337](#issue-337), [#352](#issue-352), [#357](#issue-357), [#407](#issue-407) |
| **Matching SNP frequencies and count files** | Frequency files and count sidecars belonging to the same reference release and annotation model; mixing Phase 1 and Phase 3 inputs is a concrete recurring compatibility problem. | [#46](#issue-46), [#51](#issue-51), [#61](#issue-61), [#106](#issue-106) |
| **Gene coordinate catalogs** | An available replacement for `ENSG_coord.txt`, including catalogs for different genome builds and a clear way to obtain the required gene identifiers and coordinates from a source such as GENCODE. | [#428](#issue-428), [#462](#issue-462) |
| **Working tutorial datasets** | Accessible bipolar/schizophrenia GWAS summary statistics and demonstration reference panels, with access terms and archive contents explained; a chromosome-22 example must be distinguished from a genome-wide panel. | [#179](#issue-179), [#273](#issue-273), [#282](#issue-282), [#338](#issue-338), [#403](#issue-403) |

Across these families, users also request reference sample IDs, relatedness exclusions, genome-build declarations, SNP-selection/QC provenance and annotation-definition manifests ([#33](#issue-33), [#63](#issue-63), [#224](#issue-224), [#495](#issue-495), [#505](#issue-505)).
An ancestry label alone does not resolve subpopulation matching, and an hg38 archive alone does not establish a complete compatible analysis.
LD Hub access requests concern a hosted service rather than a missing downloadable artifact ([#311](#issue-311), [#325](#issue-325), [#465](#issue-465)).

## Resource-wiki preparation

Reference resources should be selected by scientific role and compatibility, not by copying a historical filename into every workflow.
Native LDSC3 LD-score generation already prepares reference metadata, counts, overlap and regression weights together; users do not always need to acquire a separate legacy weight or frequency archive.
Imported legacy suites have a stricter, different conversion contract.

| Role | When needed | Compatibility that must be stated | Candidate source or current support | Remaining checks |
| --- | --- | --- | --- | --- |
| GWAS summary statistics | All regression analyses | Population, phenotype, allele direction, valid per-statistic N, full appropriate marker coverage and prior correction | User's GWAS; LDSC3 munging and direct legacy input | Avoid selected significant hits and unrecorded sample-size conventions. |
| PLINK reference panel or R² matrix | Computing LD scores | Ancestry, genome build, SNP identities/alleles, reference sample selection and LD window | Build current R² matrices from a suitable panel; EUR/EAS panel archives are listed in the [Gazal record](https://zenodo.org/records/10515792). | Reference sample IDs and archive contents were not independently validated. |
| Regression SNP selection | Munging and regression-row selection | SNP identity/build, alleles, coverage and exclusions | Bundled HM3 map, or an explicitly selected alternative | The map has its own attribution and coordinate-validation limits. |
| Ordinary LD scores and weights | Reusing a standard heritability/correlation reference | One compatible reference/weight pair, appropriate ancestry and complete chromosome/count families | [EUR ordinary-score mirror](https://zenodo.org/records/8182036); native generation is an alternative. | Mirror metadata identifies its origin; no numerical identity comparison was run. |
| Baseline annotations/model | Joint partitioning and query conditioning | Model version, genomic content, SNP grid, population/build, matching counts and overlap | [Baseline version history](https://zenodo.org/records/10515792/files/readme_baseline_versions.txt?download=1) and archived annotation/LD-score suites | Version names alone do not guarantee the same annotations or estimand. |
| Query annotations or gene lists | Additional query/cell-type tests | The same annotation SNP grid and reference preparation as the baseline; appropriate controls | User BED/numeric annotations or gene lists; legacy cell-type archives may require authenticated access | Legacy query LD-score suites are not generally convertible into a current query directory. |
| Gene coordinate catalog and control gene list | Gene-list queries | Declared genome build, one-based coordinates, unambiguous identifiers, explicit padding and the correct gene-expression universe | User-curated catalog in LDSC3's documented schema | No single gene catalog/control list is scientifically correct for every experiment. |
| Genetic map | cM LD windows when input CM is absent or unsuitable | Genome build, coordinate basis, units and chromosome coverage | User-supplied matching map through current map-source flags | A kb window is a scientifically different choice, not a silent substitute for missing cM. |
| Counts, frequencies and annotation overlap | Partitioned summaries; legacy conversion | Broad reference SNPs versus filtered regression SNPs; common-MAF rule; aligned annotation columns | Native generation records these; baseline legacy conversion requires a frequency suite and full annotations. | Do not calculate denominators from the printed HM3 regression rows alone. |
| Liftover chains and region exclusions | Only workflows needing these transformations | Direction, genome build, coordinate convention and region definition | Documented chain input; bundled MHC/centromere presets | Region presets and external chains have independent provenance/validation requirements. |
| Alternative-ancestry LD scores | Analyses for another ancestry | The actual target population and the method's assumptions, not merely the archive label | [Five-ancestry resource record](https://zenodo.org/records/8096589) provides a lead. | Its listing is not validation of ordinary LDSC on admixed GWAS or of LDSC3 conversion compatibility. |

The Gazal archive provides multiple versions and both EUR and EAS resources, as well as a GRCh38 archive.
The version-history file distinguishes baseline from baseline-LD and gives different recommendations for conditional cell-type coefficient tests and marginal enrichment estimation.
That is author guidance for those analyses, not a universal rule that one baseline is always best.
The existence of `GRCh38.tgz` does not establish that every hg38 query, weight and control file requested in issue #505 is present or mutually compatible.
Those internal files, reference individuals and sample-ID equivalence remain unverified because no large archive was unpacked.

For an ordinary converted suite, reference scores and weights may be the same compatible one-column suite.
For baseline-only partitioned conversion, full annotations, frequency data and count files are required; query/thin suites and arbitrary fragments are not accepted.
Current contracts are in [legacy suite conversion](../../../reference/formats/legacy-ldscores.md) and [SNP universes](../../../reference/ldscore-snp-universes.md).
Do not recommend renaming unrelated files merely to satisfy a filename check.

## Link-check interpretation

All 57 checked reference-resource URLs and their originating issues are listed below.
They were extracted from substantive issue bodies/comments using reference-resource domains and names, not from every unrelated publication, image or personal homepage mentioned in a thread.
Only headers were requested; large archives were not downloaded.

- **404:** the exact URL was missing at the check date; its scientific resource may exist elsewhere.
- **200:** the URL responded, subject to the listed caveat; this does not certify archive contents or scientific compatibility.
- **Sign-in redirect:** a cloud-console URL responded with a login page, not the reference object.
- **403:** access was forbidden from this client; not sufficient evidence that the resource was deleted.
- **400:** the exact request was rejected; a requester/billing requirement or request syntax may matter, and deletion is not established.
- **Unverified:** timeout, TLS/hostname failure, DNS failure or connection reset prevented a reliable availability judgment.

The old Broad LDSCORE landing URLs and many direct archive URLs returned 404.
The versioned Zenodo records and the separately attributed EUR mirror responded.
Historical HTTP-only pages, cloud-console links and failed network requests are kept separate from confirmed 404s.
Do not publish a mirror as authoritative solely because it is reachable; retain the original owner, archive version and role-specific compatibility checks.

| Resource URL | Observed result | Originating issues |
| --- | --- | --- |
| [https://alkesgroup.broadinstitute.org/LDSCORE/](<https://alkesgroup.broadinstitute.org/LDSCORE/>) | Missing exact URL (404) | [#321](https://github.com/bulik/ldsc/issues/321), [#370](https://github.com/bulik/ldsc/issues/370), [#382](https://github.com/bulik/ldsc/issues/382), [#384](https://github.com/bulik/ldsc/issues/384), [#386](https://github.com/bulik/ldsc/issues/386), [#392](https://github.com/bulik/ldsc/issues/392) |
| [https://alkesgroup.broadinstitute.org/LDSCORE/eur_w_ld_chr.tar.bz2/](<https://alkesgroup.broadinstitute.org/LDSCORE/eur_w_ld_chr.tar.bz2/>) | Missing exact URL (404) | [#429](https://github.com/bulik/ldsc/issues/429) |
| [https://alkesgroup.broadinstitute.org/LDSCORE/w_hm3.snplist.bz2](<https://alkesgroup.broadinstitute.org/LDSCORE/w_hm3.snplist.bz2>) | Missing exact URL (404) | [#235](https://github.com/bulik/ldsc/issues/235) |
| [https://data.broadinstitute.org/alkesgroup/LDSCORE/](<https://data.broadinstitute.org/alkesgroup/LDSCORE/>) | Missing exact URL (404) | [#26](https://github.com/bulik/ldsc/issues/26), [#43](https://github.com/bulik/ldsc/issues/43), [#46](https://github.com/bulik/ldsc/issues/46), [#51](https://github.com/bulik/ldsc/issues/51), [#61](https://github.com/bulik/ldsc/issues/61), [#63](https://github.com/bulik/ldsc/issues/63), [#103](https://github.com/bulik/ldsc/issues/103), [#106](https://github.com/bulik/ldsc/issues/106), [#235](https://github.com/bulik/ldsc/issues/235), [#370](https://github.com/bulik/ldsc/issues/370), [#377](https://github.com/bulik/ldsc/issues/377) |
| [https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase1_baseline_ldscores.tgz](<https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase1_baseline_ldscores.tgz>) | Missing exact URL (404) | [#108](https://github.com/bulik/ldsc/issues/108) |
| [https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase1_cell_type_groups.tgz](<https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase1_cell_type_groups.tgz>) | Missing exact URL (404) | [#108](https://github.com/bulik/ldsc/issues/108) |
| [https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase1_frq.tgz](<https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase1_frq.tgz>) | Missing exact URL (404) | [#108](https://github.com/bulik/ldsc/issues/108) |
| [https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baselineLD_v1.1_ldscores.tgz](<https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baselineLD_v1.1_ldscores.tgz>) | Missing exact URL (404) | [#105](https://github.com/bulik/ldsc/issues/105), [#108](https://github.com/bulik/ldsc/issues/108) |
| [https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baseline_v1.1_ldscores.tgz](<https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baseline_v1.1_ldscores.tgz>) | Missing exact URL (404) | [#96](https://github.com/bulik/ldsc/issues/96) |
| [https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_cell_type_groups.tgz](<https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_cell_type_groups.tgz>) | Missing exact URL (404) | [#108](https://github.com/bulik/ldsc/issues/108) |
| [https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_frq.tgz](<https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_frq.tgz>) | Missing exact URL (404) | [#108](https://github.com/bulik/ldsc/issues/108) |
| [https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz](<https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz>) | Missing exact URL (404) | [#190](https://github.com/bulik/ldsc/issues/190), [#414](https://github.com/bulik/ldsc/issues/414) |
| [https://data.broadinstitute.org/alkesgroup/LDSCORE/1kg_eur.tar.bz2](<https://data.broadinstitute.org/alkesgroup/LDSCORE/1kg_eur.tar.bz2>) | Missing exact URL (404) | [#403](https://github.com/bulik/ldsc/issues/403) |
| [https://data.broadinstitute.org/alkesgroup/LDSCORE/LDSC_SEG_ldscores/](<https://data.broadinstitute.org/alkesgroup/LDSCORE/LDSC_SEG_ldscores/>) | Missing exact URL (404) | [#412](https://github.com/bulik/ldsc/issues/412) |
| [https://data.broadinstitute.org/alkesgroup/LDSCORE/LDSC_SEG_ldscores/${cts_name}_1000Gv3_ldscores.tgz](<https://data.broadinstitute.org/alkesgroup/LDSCORE/LDSC_SEG_ldscores/${cts_name}_1000Gv3_ldscores.tgz>) | Missing exact URL (404) | [#450](https://github.com/bulik/ldsc/issues/450) |
| [https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2](<https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2>) | Missing exact URL (404) | [#370](https://github.com/bulik/ldsc/issues/370), [#371](https://github.com/bulik/ldsc/issues/371), [#429](https://github.com/bulik/ldsc/issues/429) |
| [https://data.broadinstitute.org/alkesgroup/LDSCORE/hapmap3_snps.tgz](<https://data.broadinstitute.org/alkesgroup/LDSCORE/hapmap3_snps.tgz>) | Missing exact URL (404) | [#96](https://github.com/bulik/ldsc/issues/96), [#322](https://github.com/bulik/ldsc/issues/322), [#331](https://github.com/bulik/ldsc/issues/331), [#337](https://github.com/bulik/ldsc/issues/337) |
| [https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz](<https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz>) | Missing exact URL (404) | [#145](https://github.com/bulik/ldsc/issues/145) |
| [https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2](<https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2>) | Missing exact URL (404) | [#45](https://github.com/bulik/ldsc/issues/45), [#96](https://github.com/bulik/ldsc/issues/96), [#235](https://github.com/bulik/ldsc/issues/235), [#370](https://github.com/bulik/ldsc/issues/370), [#371](https://github.com/bulik/ldsc/issues/371) |
| [https://data.broadinstitute.org/alkesgroup/LDSCORE/weights_hm3_no_hla.tgz](<https://data.broadinstitute.org/alkesgroup/LDSCORE/weights_hm3_no_hla.tgz>) | Missing exact URL (404) | [#108](https://github.com/bulik/ldsc/issues/108) |
| [https://data.broadinstitute.org/alkesgroup/sumstats_formatted/](<https://data.broadinstitute.org/alkesgroup/sumstats_formatted/>) | Missing exact URL (404) | [#108](https://github.com/bulik/ldsc/issues/108) |
| [https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/baselineLD_v2.2_bedfiles.tgz](<https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/baselineLD_v2.2_bedfiles.tgz>) | Missing exact URL (404) | [#233](https://github.com/bulik/ldsc/issues/233) |
| [https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/list.txt](<https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/list.txt>) | Missing exact URL (404) | [#233](https://github.com/bulik/ldsc/issues/233) |
| [http://www.broadinstitute.org/collaboration/giant/images/b/b7/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.txt.gz](<http://www.broadinstitute.org/collaboration/giant/images/b/b7/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.txt.gz>) | Forbidden (403) | [#26](https://github.com/bulik/ldsc/issues/26) |
| [http://www.broadinstitute.org/~bulik/eur_ldscores/](<http://www.broadinstitute.org/~bulik/eur_ldscores/>) | Forbidden (403) | [#23](https://github.com/bulik/ldsc/issues/23) |
| [https://www.broadinstitute.org/files/shared/mpg/hapmap3/hapmap3_r1_b36_fwd.qc.poly.tar.bz2](<https://www.broadinstitute.org/files/shared/mpg/hapmap3/hapmap3_r1_b36_fwd.qc.poly.tar.bz2>) | Forbidden (403) | [#407](https://github.com/bulik/ldsc/issues/407) |
| [https://www.broadinstitute.org/files/shared/mpg/hapmap3/hapmap3_r1_b36_fwd_consensus.qc.poly.recode.map.bz2](<https://www.broadinstitute.org/files/shared/mpg/hapmap3/hapmap3_r1_b36_fwd_consensus.qc.poly.recode.map.bz2>) | Forbidden (403) | [#407](https://github.com/bulik/ldsc/issues/407) |
| [https://www.broadinstitute.org/medical-and-population-genetics/hapmap-3](<https://www.broadinstitute.org/medical-and-population-genetics/hapmap-3>) | Forbidden (403) | [#407](https://github.com/bulik/ldsc/issues/407) |
| [https://storage.googleapis.com/broad-alkesgroup-public-requester-pays/LDSCORE/1kg_eur.tar.bz2](<https://storage.googleapis.com/broad-alkesgroup-public-requester-pays/LDSCORE/1kg_eur.tar.bz2>) | Request rejected (400) | [#321](https://github.com/bulik/ldsc/issues/321), [#402](https://github.com/bulik/ldsc/issues/402), [#403](https://github.com/bulik/ldsc/issues/403) |
| [http://ldsc.broadinstitute.org/](<http://ldsc.broadinstitute.org/>) | Unverified: timed out | [#260](https://github.com/bulik/ldsc/issues/260), [#311](https://github.com/bulik/ldsc/issues/311), [#465](https://github.com/bulik/ldsc/issues/465) |
| [http://ldsc.broadinstitute.org/about/](<http://ldsc.broadinstitute.org/about/>) | Unverified: timed out | [#57](https://github.com/bulik/ldsc/issues/57) |
| [https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baselineLD_ldscores](<https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baselineLD_ldscores>) | Unverified: [Errno 8] nodename nor servname provided, or not known | [#96](https://github.com/bulik/ldsc/issues/96) |
| [https://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/](<https://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/>) | Unverified: [SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed: Hostname mismatch, certificate is not valid for 'hgdownload.cse.ucsc.edu'. (_ssl.c:1032) | [#487](https://github.com/bulik/ldsc/issues/487) |
| [https://pan-dev.ukbb.broadinstitute.org/downloads](<https://pan-dev.ukbb.broadinstitute.org/downloads>) | Unverified: [SSL: UNEXPECTED_EOF_WHILE_READING] EOF occurred in violation of protocol (_ssl.c:1032) | [#322](https://github.com/bulik/ldsc/issues/322) |
| [http://www.1000genomes.org/data](<http://www.1000genomes.org/data>) | Responds (HEAD 200) | [#43](https://github.com/bulik/ldsc/issues/43) |
| [https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/baselineLF_v2.2.UKB.tar.gz](<https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/baselineLF_v2.2.UKB.tar.gz>) | Responds (HEAD 200) | [#384](https://github.com/bulik/ldsc/issues/384) |
| [https://console.cloud.google.com/storage/browser/_details/broad-alkesgroup-public-requester-pays/LDSCORE/1000G_Phase3_plinkfiles.tgz](<https://console.cloud.google.com/storage/browser/_details/broad-alkesgroup-public-requester-pays/LDSCORE/1000G_Phase3_plinkfiles.tgz?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22>) | Sign-in redirect (HEAD 200) | [#495](https://github.com/bulik/ldsc/issues/495) |
| [https://console.cloud.google.com/storage/browser/_details/broad-alkesgroup-public-requester-pays/LDSCORE/LDSC_SEG_ldscores/Cahoy_1000Gv3_ldscores.tgz](<https://console.cloud.google.com/storage/browser/_details/broad-alkesgroup-public-requester-pays/LDSCORE/LDSC_SEG_ldscores/Cahoy_1000Gv3_ldscores.tgz>) | Sign-in redirect (HEAD 200) | [#378](https://github.com/bulik/ldsc/issues/378) |
| [https://console.cloud.google.com/storage/browser/_details/broad-alkesgroup-public-requester-pays/LDSCORE/hapmap3_snps.tgz;tab=live_object](<https://console.cloud.google.com/storage/browser/_details/broad-alkesgroup-public-requester-pays/LDSCORE/hapmap3_snps.tgz;tab=live_object>) | Sign-in redirect (HEAD 200) | [#386](https://github.com/bulik/ldsc/issues/386) |
| [https://console.cloud.google.com/storage/browser/broad-alkesgroup-public-requester-pays](<https://console.cloud.google.com/storage/browser/broad-alkesgroup-public-requester-pays>) | Sign-in redirect (HEAD 200) | [#382](https://github.com/bulik/ldsc/issues/382), [#384](https://github.com/bulik/ldsc/issues/384), [#403](https://github.com/bulik/ldsc/issues/403), [#429](https://github.com/bulik/ldsc/issues/429), [#450](https://github.com/bulik/ldsc/issues/450) |
| [https://console.cloud.google.com/storage/browser/broad-alkesgroup-public-requester-pays/LDSCORE/LDSC_SEG_ldscores](<https://console.cloud.google.com/storage/browser/broad-alkesgroup-public-requester-pays/LDSCORE/LDSC_SEG_ldscores>) | Sign-in redirect (HEAD 200) | [#412](https://github.com/bulik/ldsc/issues/412) |
| [https://console.cloud.google.com/storage/browser/broad-alkesgroup-public-requester-pays/LDSCORE](<https://console.cloud.google.com/storage/browser/broad-alkesgroup-public-requester-pays/LDSCORE?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22>) | Sign-in redirect (HEAD 200) | [#235](https://github.com/bulik/ldsc/issues/235), [#470](https://github.com/bulik/ldsc/issues/470) |
| [https://doi.org/10.5281/zenodo.7768714](<https://doi.org/10.5281/zenodo.7768714>) | Responds (HEAD 200) | [#384](https://github.com/bulik/ldsc/issues/384) |
| [https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/README_new_data_location.txt](<https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/README_new_data_location.txt>) | Responds (HEAD 200) | [#382](https://github.com/bulik/ldsc/issues/382) |
| [https://zenodo.org/record/7768714](<https://zenodo.org/record/7768714>) | Responds (HEAD 200) | [#389](https://github.com/bulik/ldsc/issues/389) |
| [https://zenodo.org/records/10515792](<https://zenodo.org/records/10515792>) | Responds (HEAD 200) | [#429](https://github.com/bulik/ldsc/issues/429), [#437](https://github.com/bulik/ldsc/issues/437), [#476](https://github.com/bulik/ldsc/issues/476), [#479](https://github.com/bulik/ldsc/issues/479), [#499](https://github.com/bulik/ldsc/issues/499), [#505](https://github.com/bulik/ldsc/issues/505) |
| [https://zenodo.org/records/10515792/files/1000G_Phase3_EAS_weights_hm3_no_MHC.tgz](<https://zenodo.org/records/10515792/files/1000G_Phase3_EAS_weights_hm3_no_MHC.tgz?download=1>) | Responds (HEAD 200) | [#376](https://github.com/bulik/ldsc/issues/376) |
| [https://zenodo.org/records/10515792/files/1000G_Phase3_ldscores.tgz](<https://zenodo.org/records/10515792/files/1000G_Phase3_ldscores.tgz?download=1>) | Responds (HEAD 200) | [#429](https://github.com/bulik/ldsc/issues/429) |
| [https://zenodo.org/records/10515792/files/GRCh38.tgz](<https://zenodo.org/records/10515792/files/GRCh38.tgz>) | Responds (HEAD 200) | [#264](https://github.com/bulik/ldsc/issues/264) |
| [https://zenodo.org/records/10515792/files/GRCh38.tgz](<https://zenodo.org/records/10515792/files/GRCh38.tgz?download=1>) | Responds (HEAD 200) | [#470](https://github.com/bulik/ldsc/issues/470) |
| [https://zenodo.org/records/7768714](<https://zenodo.org/records/7768714>) | Responds (HEAD 200) | [#414](https://github.com/bulik/ldsc/issues/414), [#474](https://github.com/bulik/ldsc/issues/474) |
| [https://zenodo.org/records/7768714/files/1000G_Phase3_EAS_weights_hm3_no_MHC.tgz](<https://zenodo.org/records/7768714/files/1000G_Phase3_EAS_weights_hm3_no_MHC.tgz?download=1>) | Responds (HEAD 200) | [#425](https://github.com/bulik/ldsc/issues/425), [#470](https://github.com/bulik/ldsc/issues/470) |
| [https://zenodo.org/records/7768714/files/1000G_Phase3_weights_hm3_no_MHC.tgz](<https://zenodo.org/records/7768714/files/1000G_Phase3_weights_hm3_no_MHC.tgz?download=1>) | Responds (HEAD 200) | [#425](https://github.com/bulik/ldsc/issues/425) |
| [https://zenodo.org/records/7768714/files/baselineLD_v2.2_bedfiles.tgz](<https://zenodo.org/records/7768714/files/baselineLD_v2.2_bedfiles.tgz?download=1>) | Responds (HEAD 200) | [#470](https://github.com/bulik/ldsc/issues/470) |
| [https://zenodo.org/records/8096589](<https://zenodo.org/records/8096589>) | Responds (HEAD 200) | [#503](https://github.com/bulik/ldsc/issues/503) |
| [https://zenodo.org/records/8182036](<https://zenodo.org/records/8182036>) | Responds (HEAD 200) | [#421](https://github.com/bulik/ldsc/issues/421), [#482](https://github.com/bulik/ldsc/issues/482) |
| [https://zenodo.org/records/8292725](<https://zenodo.org/records/8292725>) | Responds (HEAD 200) | [#233](https://github.com/bulik/ldsc/issues/233) |

## Recommendations by kind

- **Code:** keep conversion, genome-build and SNP-grid checks; resource availability cannot be repaired by silently relaxing them.
- **Documentation:** curate a short matrix of analysis goal, resource role, ancestry, genome build and baseline/control pairing, with verified links and dates.
- **Features:** a new automatic downloader or resource database is unnecessary at this stage.
- **Resources:** prefer stable owner/DOI landing pages; record mirrors and unresolved licenses separately, and validate archive contents before claiming a ready-to-run combination.

## Category summary

54 issues and 0 FAQ entries.

| Severity | Source | Concern | LDSC3 disposition | Evidence level |
| --- | --- | --- | --- | --- |
| High | [#63](#issue-63) | Tutorial reference genome build is unstated | Partial coverage | Source inspection |
| High | [#322](#issue-322) | HM3 provenance and rsID evolution | Partial coverage | Source inspection |
| High | [#357](#issue-357) | hg38 HM3 map | Partial coverage | Source inspection |
| High | [#407](#issue-407) | AFR HM3 and weight roles | Partial coverage | Source inspection |
| High | [#470](#issue-470) | hg38 resource compatibility | Partial coverage | Source inspection |
| High | [#51](#issue-51) | Mixed reference releases mistaken for memory failure | Resource-only concern | Source inspection |
| High | [#235](#issue-235) | HM3/reference hosting migrations | Resource-only concern | Source inspection |
| High | [#248](#issue-248) | East Asian baseline/query/weights pairing | Resource-only concern | Source inspection |
| High | [#337](#issue-337) | Missing EAS HM3 regression-SNP list | Resource-only concern | Source inspection |
| High | [#352](#issue-352) | EAS SNP-grid compatibility | Resource-only concern | Source inspection |
| High | [#376](#issue-376) | Ancestry/species-specific weights | Resource-only concern | Source inspection |
| High | [#479](#issue-479) | CHB versus EAS reference | Resource-only concern | Source inspection |
| High | [#495](#issue-495) | VCF-to-PLINK QC provenance | Resource-only concern | Source inspection |
| High | [#503](#issue-503) | Additional ancestry and multi-ancestry resources | Resource-only concern | Source inspection |
| High | [#505](#issue-505) | Complete EUR hg38 resource tuple | Resource-only concern | Source inspection |
| High | [#264](#issue-264) | GRCh38 summary statistics with GRCh37 resources | Apparently prevented | Source inspection |
| High | [#331](#issue-331) | Genome build of downloaded PLINK/HM3 files | Apparently prevented | Source inspection |
| Moderate | [#105](#issue-105) | Historical annotation header mixes spaces and tabs | Partial coverage | Source inspection |
| Moderate | [#384](#issue-384) | HM3/archive discovery | Partial coverage | Source inspection |
| Moderate | [#386](#issue-386) | HM3 archive location | Partial coverage | Source inspection |
| Moderate | [#428](#issue-428) | Gene-coordinate catalog provenance | Partial coverage | Source inspection |
| Moderate | [#462](#issue-462) | Gene-coordinate source | Partial coverage | Source inspection |
| Moderate | [#33](#issue-33) | Reference-panel sample membership | Resource-only concern | Source inspection |
| Moderate | [#45](#issue-45) | Relocated HapMap3 restriction list | Resource-only concern | Source inspection |
| Moderate | [#46](#issue-46) | Files hidden inside resource archives | Resource-only concern | Source inspection |
| Moderate | [#61](#issue-61) | Unclear archive contents in partitioning tutorial | Resource-only concern | Source inspection |
| Moderate | [#106](#issue-106) | Cell-type examples hidden under changed archive names | Resource-only concern | Source inspection |
| Moderate | [#169](#issue-169) | Local package and LD Hub use different SNP resources | Resource-only concern | Source inspection |
| Moderate | [#179](#issue-179) | PGC tutorial download links | Resource-only concern | Source inspection |
| Moderate | [#201](#issue-201) | GTEx brain versus multi-tissue annotations | Resource-only concern | Source inspection |
| Moderate | [#208](#issue-208) | Broken genetic-correlation tutorial | Resource-only concern | Source inspection |
| Moderate | [#224](#issue-224) | Baseline annotation acronym manifest | Resource-only concern | Source inspection |
| Moderate | [#231](#issue-231) | Repeated PGC direct-download404 | Resource-only concern | Source inspection |
| Moderate | [#273](#issue-273) | PGC tutorial requires access form | Resource-only concern | Source inspection |
| Moderate | [#282](#issue-282) | Tutorial dataset link and mirror fail | Resource-only concern | Source inspection |
| Moderate | [#311](#issue-311) | LD Hub service unavailable | Resource-only concern | Source inspection |
| Moderate | [#325](#issue-325) | LD Hub outage and alternatives | Resource-only concern | Source inspection |
| Moderate | [#330](#issue-330) | Redirected PGC URL still404 | Resource-only concern | Source inspection |
| Moderate | [#338](#issue-338) | Duplicate PGC tutorial404 | Resource-only concern | Source inspection |
| Moderate | [#370](#issue-370) | Intermittent Broad resource access | Resource-only concern | Source inspection |
| Moderate | [#371](#issue-371) | Moved EUR/HM3 resources | Resource-only concern | Source inspection |
| Moderate | [#377](#issue-377) | EUR/HM3 unavailable | Resource-only concern | Source inspection |
| Moderate | [#378](#issue-378) | Requester-pays billing | Resource-only concern | Source inspection |
| Moderate | [#382](#issue-382) | EAS archive sharing | Resource-only concern | Source inspection |
| Moderate | [#402](#issue-402) | EAS resource access | Resource-only concern | Source inspection |
| Moderate | [#403](#issue-403) | Tutorial panel versus genome-wide panel | Resource-only concern | Source inspection |
| Moderate | [#409](#issue-409) | Requester-pays availability | Resource-only concern | Source inspection |
| Moderate | [#412](#issue-412) | Moved cell-type archives | Resource-only concern | Source inspection |
| Moderate | [#414](#issue-414) | EUR PLINK panel location | Resource-only concern | Source inspection |
| Moderate | [#429](#issue-429) | Replacement archive ambiguity | Resource-only concern | Source inspection |
| Moderate | [#450](#issue-450) | Cahoy tutorial archive moved | Resource-only concern | Source inspection |
| Moderate | [#465](#issue-465) | LD Hub service unavailable | Resource-only concern | Source inspection |
| Moderate | [#499](#issue-499) | EAS archive role versus basename | Resource-only concern | Source inspection |
| Low | [#219](#issue-219) | Canonical paper URL | Resource-only concern | Source inspection |

## Individual assessments

<a id="issue-63"></a>

### #63: Genome build of tutorial data?

[Original issue](https://github.com/bulik/ldsc/issues/63) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **2**.

Workflow tags: `ldscore`.

**Trigger and harm.** Users infer hg19 from a few coordinates because release metadata omit it.
Wrong-build interval projection or LD windows can yield scientifically incorrect annotations.

**Legacy disposition.** Maintainer confirms historical resources use hg19 and notes old rsID joins; that does not validate coordinate mixing.

**LDSC3 assessment — Partial coverage.** LDSC3 adds explicit hg19/hg38 policy, inference and coordinate-family identity checks.
Resource release declarations remain necessary; rsID-family matching can join without proving that BED intervals, gene coordinates or excluded regions have the same genome build.
The old comment that absolute coordinates matter little must not be reused for annotation workflows.

**Recommended prevention.** Require genome build and coordinate basis in every resource entry, and explain that rsID joins do not authorize mismatched interval coordinates.

**Evidence and verification.** source-inspected
[src/ldsc/_genome/build_policy.py](<../../../../src/ldsc/_genome/build_policy.py>); [docs/reference/snp-identity-and-genome-build.md](<../../../reference/snp-identity-and-genome-build.md>); [docs/reference/formats/bed-intervals.md](<../../../reference/formats/bed-intervals.md>).

<a id="issue-322"></a>

### #322: HapMap3 variants

[Original issue](https://github.com/bulik/ldsc/issues/322) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **3**.

Workflow tags: `munge-sumstats`, `ldscore`.

**Trigger and harm.** Liftover of an old hg18 HM3 release loses variants retained by other tools.
Assuming every HM3 list is identical can create missing/wrong matches across biobanks.

**Legacy disposition.** No maintainer provenance supplied; users report their own liftover and external leads without validation.

**LDSC3 assessment — Partial coverage.** LDSC3 bundles an allele-aware HM3 map with hg19/hg38 positions and separate build-inference references.
This supplies an operational consistent default but does not retrospectively establish the exact curation of every legacy HM3 list; attribution notes explicitly retain unresolved provenance/reuse questions.

**Recommended prevention.** Publish known source/release/filtering facts, distinguish renamed rsIDs from unmappable loci and state unresolved provenance rather than inventing a lineage.

**Evidence and verification.** source-inspected
[src/ldsc/data/ATTRIBUTION.txt](<../../../../src/ldsc/data/ATTRIBUTION.txt>); [src/ldsc/data/readme.txt](<../../../../src/ldsc/data/readme.txt>); [src/ldsc/_genome/hm3.py](<../../../../src/ldsc/_genome/hm3.py>).

<a id="issue-357"></a>

### #357: w_hm3.snplist for GRCh38

[Original issue](https://github.com/bulik/ldsc/issues/357) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **3**.
The API reports 4 comments; one was unavailable through both retrieval methods.

Workflow tags: `munge-sumstats`, `ldscore`.

**Trigger and harm.** Requests a GRCh38 allele-merge SNP list; comments construct one by rsID and PLINK merging. A reused hg19 list or unsafe rsID join can silently mismatch alleles or coordinates.

**Legacy disposition.** Community workaround only; 3 comments were accessible although GitHub reports 4, so the thread is incompletely observable.

**LDSC3 assessment — Partial coverage.** LDSC3 bundles HM3 positions for both supported genome builds and offers allele-aware matching.
The map provenance and revalidation limits still matter, and the thread recipe is not a verified substitute for a curated resource.

**Recommended prevention.** Document the packaged map role, genome-build selection and validation limits; avoid requiring users to rebuild an allele list by unqualified rsID matching.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_genome/hm3.py](<../../../../src/ldsc/_genome/hm3.py>); [src/ldsc/data/ATTRIBUTION.txt](<../../../../src/ldsc/data/ATTRIBUTION.txt>); [docs/dev/validation/hm3-map-coordinates.md](<../hm3-map-coordinates.md>).

<a id="issue-407"></a>

### #407: Suitable w_hm3.snplist for AFR?

[Original issue](https://github.com/bulik/ldsc/issues/407) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **10**.

Workflow tags: `munge-sumstats`, `ldscore`, `rg`.

**Trigger and harm.** Requests an AFR matching list; comments disagree about ancestry specificity and confuse LD contributors with weight SNPs.
Using a European LD panel or restricting LD contributors to HM3 can change the scientific analysis.

**Legacy disposition.** Comments provide historical HapMap files, a build-liftover requirement and a linked maintainer statement; no single validated AFR bundle is established.

**LDSC3 assessment — Partial coverage.** LDSC3 separates reference SNPs from regression SNPs and computes weights accordingly. A generic HM3 identity list does not supply ancestry-matched LD or prove portability to admixed data.

**Recommended prevention.** Explain the two SNP sets with a diagram/table, require ancestry-matched LD, and qualify historical b36 files and community recipes before resource-wiki inclusion.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [src/ldsc/_genome/hm3.py](<../../../../src/ldsc/_genome/hm3.py>); [src/ldsc/_ldscore/calculator.py](<../../../../src/ldsc/_ldscore/calculator.py>).

<a id="issue-470"></a>

### #470: LDSC on build38 data

[Original issue](https://github.com/bulik/ldsc/issues/470) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **5**.

Workflow tags: `munge-sumstats`, `ldscore`, `rg`.

**Trigger and harm.** Requests build38 analysis without lifting GWAS to build37; later comments link GRCh38.tgz. rsID-only joins can conceal incompatible coordinate-based annotations, maps or allele identities.

**Legacy disposition.** Comments first describe private rebuilt resources, then point to a public archive; no complete compatible tuple is verified in-thread.

**LDSC3 assessment — Partial coverage.** LDSC3 supports hg19/hg38, inference and explicit conversions, but every positional resource must agree. A public archive label does not establish which baseline, weights, gene catalog or panel should be paired.

**Recommended prevention.** Publish a verified hg38 resource matrix and explain when rsID matching is insufficient; qualify the ancestry-specific example resources rather than copying them universally.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/reference/snp-identity-and-genome-build.md](<../../../reference/snp-identity-and-genome-build.md>); [src/ldsc/_genome/build_inference.py](<../../../../src/ldsc/_genome/build_inference.py>); [src/ldsc/_ldscore/legacy_converter.py](<../../../../src/ldsc/_ldscore/legacy_converter.py>).

<a id="issue-51"></a>

### #51: link not found

[Original issue](https://github.com/bulik/ldsc/issues/51) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **7**.

Workflow tags: `baseline-partitioned-h2`.

**Trigger and harm.** A temporary resource outage is followed by Phase 1 frequency files combined with Phase 3 LD scores. A generic annotation error can prompt wasted memory tuning while resource identities are inconsistent.

**Legacy disposition.** Maintainer identifies release mismatch and uploads matching Phase 3 frequency files; outage itself was temporary.

**LDSC3 assessment — Resource-only concern.** LDSC3 validates saved artifact schemas/counts and native generation derives MAF from the selected panel.
It cannot prove ancestry/release coherence from arbitrary filenames alone, particularly after legacy conversion. A successful URL response also does not establish that all supplied files belong together.

**Recommended prevention.** Make resource bundles list release, ancestry, build and companion requirements; diagnose identity/schema before suggesting more RAM.

**Evidence and verification.** source-inspected
[src/ldsc/_ldscore/legacy_converter.py](<../../../../src/ldsc/_ldscore/legacy_converter.py>); [docs/reference/formats/legacy-ldscores.md](<../../../reference/formats/legacy-ldscores.md>); [docs/user-guide/reference-data.md](<../../../user-guide/reference-data.md>).

<a id="issue-235"></a>

### #235: HM3 SNPlist is gone or moved

[Original issue](https://github.com/bulik/ldsc/issues/235) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **7**.

Workflow tags: `munge-sumstats`, `h2`, `rg`.

**Trigger and harm.** Broad resource directory disappears and the HM3 list URL fails.
Common preprocessing and reference acquisition become unavailable across many workflows.

**Legacy disposition.** Maintainers report hosting relocation and redirects in2020; later comments report renewed link/access problems, so closure is not permanent availability proof.

**LDSC3 assessment — Resource-only concern.** LDSC3 bundles its HM3 map, reducing this exact preprocessing dependency, but external LD-score directories and reference panels still require maintained provenance and access instructions. A requester-pays URL is an access-policy change, not necessarily an absent file.

**Recommended prevention.** Inventory official landing pages and mirrors with current link/access status; identify which minimal resources are bundled.

**Evidence and verification.** source-inspected
[src/ldsc/data/hm3_curated_map.tsv.gz](<../../../../src/ldsc/data/hm3_curated_map.tsv.gz>); [src/ldsc/data/ATTRIBUTION.txt](<../../../../src/ldsc/data/ATTRIBUTION.txt>); [Source](<https://github.com/bulik/ldsc/issues/235>).

<a id="issue-248"></a>

### #248: Files chosen for cell-type specific analysis on EAS population

[Original issue](https://github.com/bulik/ldsc/issues/248) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** A user needs the EAS equivalent of the documented EUR cell-type analysis resources.
Mixing European LD with East Asian GWAS changes the model and can bias interpretation.

**Legacy disposition.** Closed with no public answer or identified fix.

**LDSC3 assessment — Resource-only concern.** LDSC3 accepts compatible panels/directories but does not automatically supply or scientifically validate an EAS resource suite.
Query gene definitions may be reusable, whereas LD scores, regression weights, frequencies and common-SNP counts must be generated consistently for the intended ancestry.

**Recommended prevention.** Catalog complete ancestry/build-matched bundles and a regeneration route; distinguish transferable biological annotations from population-specific LD.

**Evidence and verification.** source-inspected
[docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [docs/user-guide/workflows/cell-type-enrichment.md](<../../../user-guide/workflows/cell-type-enrichment.md>).

<a id="issue-337"></a>

### #337: Can't find EAS HapMap3 SNPs in alkesgroup repository

[Original issue](https://github.com/bulik/ldsc/issues/337) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **2**.

Workflow tags: `ldscore`, `ldsc-seg`.

**Trigger and harm.** Distributed EAS baseline scores use a SNP list differing from generic hapmap3_snps.tgz.
Users cannot reproduce compatible query rows or confirm the resource build.

**Legacy disposition.** Open with repeated requests and no verified published file.

**LDSC3 assessment — Resource-only concern.** LDSC3 can select common regression SNPs and generate query/baseline scores together, but cannot recover an unpublished historical list from an internal Price-lab path. A generic HM3 map is not proof of identical historical EAS filtering or MHC exclusions.

**Recommended prevention.** Publish exact available regression-SNP selection/provenance for new bundles and provide a coherent regeneration path when historical inputs are unavailable.

**Evidence and verification.** source-inspected
[docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [src/ldsc/data/readme.txt](<../../../../src/ldsc/data/readme.txt>); [Source](<https://github.com/bulik/ldsc/issues/337>).

<a id="issue-352"></a>

### #352: ValueError: LD Scores for concatenation must have identical SNP columns.

[Original issue](https://github.com/bulik/ldsc/issues/352) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `ldscore`.

**Trigger and harm.** Requests the SNP list corresponding to EAS LD-score files after a concatenation mismatch.
Combining mismatched baseline/query rows can invalidate a fitted model.

**Legacy disposition.** Open with no answer or verified resource pairing.

**LDSC3 assessment — Resource-only concern.** LDSC3 can read ancestry-matched inputs, but its packaged HM3 restriction does not prove that an EAS LD-score suite shares a compatible annotation SNP grid and weights.
Conversion validates a suite rather than supplying its missing scientific provenance.

**Recommended prevention.** Publish complete EAS resource bundles with row identifiers, ancestry, build, count conventions and compatible weights; diagnose missing versus contradictory SNPs separately.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_ldscore/legacy_converter.py](<../../../../src/ldsc/_ldscore/legacy_converter.py>); [docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>).

<a id="issue-376"></a>

### #376: weights for non-European populations

[Original issue](https://github.com/bulik/ldsc/issues/376) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **4**.

Workflow tags: `ldscore`, `h2`.

**Trigger and harm.** Requests non-European weights; comments first supply EAS weights, then reveal a chicken analysis.
EAS weights cannot validate another species and ancestry-mismatched weights can change precision/model behavior.

**Legacy disposition.** An EAS Zenodo file and weight-estimation tutorial are linked; no nonhuman validation or code fix is supplied.

**LDSC3 assessment — Resource-only concern.** LDSC3 computes regression-row weights from its reference inputs and validates imported suites, but human chromosome/build defaults make chicken analysis unsupported.
This thread spans resource availability and a separate species limitation.

**Recommended prevention.** Separate ancestry-matched human resource instructions from unsupported nonhuman requests; document the weight SNP set and LD window used.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [src/ldsc/_ldscore/legacy_converter.py](<../../../../src/ldsc/_ldscore/legacy_converter.py>); [src/ldsc/_core/chromosomes.py](<../../../../src/ldsc/_core/chromosomes.py>).

<a id="issue-479"></a>

### #479: How to make customize `panel` ?

[Original issue](https://github.com/bulik/ldsc/issues/479) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `ldscore`, `baseline-partitioned-h2`.

**Trigger and harm.** Requests CHB-specific scores/weights and asks whether baseline BEDs are European. A broad EAS panel may not exactly match CHB, and genomic annotation sources are distinct from population-derived LD/frequency features.

**Legacy disposition.** Comment points to EAS archives but does not establish equivalence to a CHB-only panel.

**LDSC3 assessment — Resource-only concern.** LDSC3 can score a chosen PLINK panel and supplied annotations, but baselineLD includes population-dependent features beyond simple BED membership.
Reusing EUR values wholesale is not justified by coordinate compatibility.

**Recommended prevention.** Separate reusable functional intervals from ancestry-dependent LD/MAF annotations and document panel composition, genetic maps, weights and counts when offering CHB/EAS resources.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_ldscore/calculator.py](<../../../../src/ldsc/_ldscore/calculator.py>); [docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [docs/reference/formats/bed-intervals.md](<../../../reference/formats/bed-intervals.md>).

<a id="issue-495"></a>

### #495: Quality Control when converting 1000 Genome Project Phase 3 vcf data into plink file

[Original issue](https://github.com/bulik/ldsc/issues/495) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `ldscore`.

**Trigger and harm.** A converted Phase3 VCF has far more variants than the distributed reference panel.
Reproducing a different variant universe changes LD, counts and downstream comparisons.

**Legacy disposition.** Open with no historical QC manifest or verified explanation.

**LDSC3 assessment — Resource-only concern.** LDSC3 validates trios and records its own reference-SNP funnel, but cannot infer which release, SNP/indel, frequency or relatedness filters produced someone else's panel.

**Recommended prevention.** Resource curation should include exact panel source and sample/variant QC, keeping raw VCF scope separate from biallelic reference SNP and regression sets.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_runtime/paths.py:454](<../../../../src/ldsc/_runtime/paths.py#L454>); [docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [tests/ldscore/test_reference_snp_funnel.py](<../../../../tests/ldscore/test_reference_snp_funnel.py>).

<a id="issue-503"></a>

### #503: LDSC for non EUR/EAS ancestries

[Original issue](https://github.com/bulik/ldsc/issues/503) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `baseline-partitioned-h2`.

**Trigger and harm.** Requests non-EUR/EAS baselines and a multi-ancestry partitioning workflow.
An available hg38 panel or ancestry mixture is not automatically a calibrated multi-ancestry LDSC model.

**Legacy disposition.** Comment links Zenodo8096589 but explicitly doubts its fit; no validated multi-ancestry remedy is supplied.

**LDSC3 assessment — Resource-only concern.** LDSC3 can compute scores from a selected panel but does not implement a dedicated multi-ancestry estimator.
Population-dependent baselineLD features, admixture and matching GWAS LD require explicit scientific validation.

**Recommended prevention.** List resource availability separately from methodological support and qualify the linked archive until ancestry/build/annotation roles are verified; do not pool panels by convenience.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [src/ldsc/_ldscore/calculator.py](<../../../../src/ldsc/_ldscore/calculator.py>); [Source](<https://zenodo.org/records/8096589>).

<a id="issue-505"></a>

### #505: Where can I find European-ancestry (EUR) hg38/GRCh38 LDSC reference files?

[Original issue](https://github.com/bulik/ldsc/issues/505) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `ldscore`, `h2`, `rg`, `ldsc-seg`.

**Trigger and harm.** Requests compatible hg38 scores, weights, baselinev1.2, Cahoy queries and panel sample identities for comparison with hg19. A bare GRCh38 archive link does not establish a complete compatible analysis or matched reference people.

**Legacy disposition.** Open without an answer; no resource pairing or sample-equivalence evidence is available in the thread.

**LDSC3 assessment — Resource-only concern.** LDSC3 supports hg38 and enforces declared scientific identities, but resource curation must still verify contents, ancestry, sample/QC provenance, baseline/query compatibility and LD window.
It cannot infer that hg19/hg38 panels used the same people.

**Recommended prevention.** Make this a resource-wiki priority: publish a verified workflow-by-build matrix and explicit unknowns; avoid claiming Cahoy hg38 or panel sample equivalence until direct evidence exists.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/reference/snp-identity-and-genome-build.md](<../../../reference/snp-identity-and-genome-build.md>); [src/ldsc/_ldscore/legacy_converter.py](<../../../../src/ldsc/_ldscore/legacy_converter.py>); [src/ldsc/data/ATTRIBUTION.txt](<../../../../src/ldsc/data/ATTRIBUTION.txt>).

<a id="issue-264"></a>

### #264: What if GWAS summary statistics are from a different reference genome build?

[Original issue](https://github.com/bulik/ldsc/issues/264) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **2**.

Workflow tags: `munge-sumstats`, `ldsc-seg`.

**Trigger and harm.** Cell-type analysis combines GWAS hg38 coordinates with hg19 references.
Mismatched coordinates lead to lost/wrong SNP matches and annotation projection.

**Legacy disposition.** A commenter supplies GRCh38.tgz on Zenodo10515792; no universal compatibility validation or code fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 resolves genome-build declarations, supports explicit hg19/hg38 liftover for munging and checks compatible artifact provenance.
It still requires the user to choose a coherent reference bundle and validate chain direction/annotation catalogs; rsID matching alone does not make all annotation resources build-independent.

**Recommended prevention.** Provide a choice between a verified native-hg38 bundle and explicit liftover into a complete hg19 workflow, with provenance and counts.

**Evidence and verification.** source-inspected
[src/ldsc/_genome/build_policy.py](<../../../../src/ldsc/_genome/build_policy.py>); [docs/reference/snp-identity-and-genome-build.md](<../../../reference/snp-identity-and-genome-build.md>).

<a id="issue-331"></a>

### #331: 1000G_Phase3_plinkfiles.tgz - GRCh37 or GRCh38

[Original issue](https://github.com/bulik/ldsc/issues/331) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `ldscore`.

**Trigger and harm.** Title names Phase3 while body asks about Phase1 and comments ask about HM3 genome build.
Filename ambiguity can combine incompatible coordinates and annotation resources.

**Legacy disposition.** Open without an authoritative answer; no code fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 separates explicit genome-build assertions, inference and artifact provenance.
The issue contradicts itself about the resource release, so the archive contents/source metadata must be checked rather than assuming the title identifies the file.
HM3 rsID membership is distinct from position-specific files.

**Recommended prevention.** Catalog exact release, build and sample ancestry; demonstrate infer-genome-build as a check and explain its unknown/conflict outcomes.

**Evidence and verification.** source-inspected
[src/ldsc/_genome/build_policy.py](<../../../../src/ldsc/_genome/build_policy.py>); [docs/reference/snp-identity-and-genome-build.md](<../../../reference/snp-identity-and-genome-build.md>).

<a id="issue-105"></a>

### #105: Mixed delimiters in header of baselineLD_v1.1 annotation

[Original issue](https://github.com/bulik/ldsc/issues/105) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **3**.

Workflow tags: `baseline-partitioned-h2`.

**Trigger and harm.** A baselineLD v1.1 header has space-separated MAF columns above tab-separated data.
External tools using strict TSV parsing can misread column counts even when legacy LDSC succeeds.

**Legacy disposition.** Maintainers retain the historical file for reproducibility and promise cleaner future releases, not an in-place repair.

**LDSC3 assessment — Partial coverage.** LDSC3 accepts documented whitespace annotation input, but generic readers and user curation still need to inspect header consistency.
Treating every whitespace-separated table as a strict TSV is incorrect.
Normalizing a downloaded copy is reasonable only with preserved source provenance and verified row/column counts, not replacing the historical original.

**Recommended prevention.** Record the known format caveat in resource entries and provide a validated local-normalization recipe without modifying published originals.

**Evidence and verification.** source-inspected
[src/ldsc/_annotation/parsing.py](<../../../../src/ldsc/_annotation/parsing.py>); [docs/reference/formats/column-schema.md](<../../../reference/formats/column-schema.md>); [Source](<https://github.com/bulik/ldsc/issues/105#issuecomment-398909510>).

<a id="issue-384"></a>

### #384: How to obtain the file "w_hm3.snplist.bz2"

[Original issue](https://github.com/bulik/ldsc/issues/384) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **7**.

Workflow tags: `munge-sumstats`, `h2`, `rg`.

**Trigger and harm.** Users confuse a compressed EUR archive, its extracted directory and a separately hosted HM3 list.
They cannot locate usable inputs and may combine unrelated baseline and weight releases.

**Legacy disposition.** Comments provide Colorado, cloud and Zenodo locations and a CELLECT HM3 copy; these are resource pointers, not code fixes.

**LDSC3 assessment — Partial coverage.** The packaged HM3 map removes the routine external-list requirement in LDSC3.
Converted legacy LD scores still require complete suites; an archive and an extracted folder can supply the same role without having the same path syntax.

**Recommended prevention.** Explain archive extraction versus file selection and publish resource-role recipes rather than making users search for one obsolete basename.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_genome/hm3.py](<../../../../src/ldsc/_genome/hm3.py>); [docs/reference/input-paths.md](<../../../reference/input-paths.md>); [docs/reference/formats/legacy-ldscores.md](<../../../reference/formats/legacy-ldscores.md>).

<a id="issue-386"></a>

### #386: Where to obtain "hapmap3_snps.tgz"?

[Original issue](https://github.com/bulik/ldsc/issues/386) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `munge-sumstats`, `ldscore`.

**Trigger and harm.** Requests hapmap3_snps.tgz after the old host cannot provide it. A filename-specific tutorial dependency blocks preparation.

**Legacy disposition.** Comment links a requester-pays object; no accessible free verified substitute is established.

**LDSC3 assessment — Partial coverage.** LDSC3 includes its default HM3 map, while users may explicitly supply another regression SNP set.
That default is not a general substitute for every historical archive member without checking its purpose.

**Recommended prevention.** Describe whether a workflow needs SNP identifiers, allele matching or build coordinates, and point to the packaged table when it satisfies that role.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_genome/hm3.py](<../../../../src/ldsc/_genome/hm3.py>); [docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>).

<a id="issue-428"></a>

### #428: where can i find ENSG_coord.txt or similiar plain file in different build?

[Original issue](https://github.com/bulik/ldsc/issues/428) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `ldscore`, `ldsc-seg`.

**Trigger and harm.** The legacy default ENSG_coord.txt is unavailable for a desired genome build.
An improvised coordinate conversion can shift intervals or map incompatible gene identifiers.

**Legacy disposition.** Comments suggest GENCODE GTF and a community conversion script; neither verifies the resulting file for this user.

**LDSC3 assessment — Partial coverage.** LDSC3 requires an explicit gene coordinate catalog with genome_build and one-based coordinates, removing the hidden default file.
It still needs a trustworthy source release and a correct GTF-to-gene reduction.

**Recommended prevention.** Provide a small versioned catalog recipe with chromosome, identifier-version, transcript/gene and coordinate conventions; validate the output instead of silently downloading a generic file.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/reference/formats/gene-lists.md](<../../../reference/formats/gene-lists.md>); [src/ldsc/_annotation/gene_lists.py](<../../../../src/ldsc/_annotation/gene_lists.py>); [docs/reference/snp-identity-and-genome-build.md](<../../../reference/snp-identity-and-genome-build.md>).

<a id="issue-462"></a>

### #462: How to generate one ENSG_coord.txt file for the make_annot.py ?

[Original issue](https://github.com/bulik/ldsc/issues/462) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `ldscore`, `ldsc-seg`.

**Trigger and harm.** Requests how to generate the legacy ENSG_coord file.
Wrong gene release or coordinate convention shifts SNP membership.

**Legacy disposition.** Open with no answer or verified recipe.

**LDSC3 assessment — Partial coverage.** LDSC3 specifies an explicit catalog containing gene identifiers, names, chromosome, start/end and genome_build.
It validates coordinates but cannot decide which annotation release or transcript aggregation is scientifically appropriate.

**Recommended prevention.** Provide a versioned GENCODE/Ensembl-to-gene catalog example with one-based coordinates, identifier-version policy and a tiny validated output.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/reference/formats/gene-lists.md](<../../../reference/formats/gene-lists.md>); [src/ldsc/_annotation/gene_lists.py](<../../../../src/ldsc/_annotation/gene_lists.py>).

<a id="issue-33"></a>

### #33: Number of EUR individuals in 1000 Genomes Project

[Original issue](https://github.com/bulik/ldsc/issues/33) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `ldscore`.

**Trigger and harm.** A historical panel has 378 EUR individuals despite 379 in a public release.
Unrecorded relatedness exclusions impede exact reproduction and reference matching.

**Legacy disposition.** Maintainer recalls excluding one member of a cousin pair; no authoritative sample manifest is provided.

**LDSC3 assessment — Resource-only concern.** LDSC3 can restrict individuals in PLINK preparation, but it cannot reconstruct historical membership from a panel label.
Modern panels may legitimately use different sample sets.
Resource documentation should identify release, ancestry, included samples and exclusions rather than copying 378/379 as a requirement for all LDSC analyses.

**Recommended prevention.** Curate sample manifests or source statements alongside resource releases and disclose when exact historical membership is unavailable.

**Evidence and verification.** source-inspected
[src/ldsc/_refpanel/plink.py](<../../../../src/ldsc/_refpanel/plink.py>); [docs/user-guide/reference-data.md](<../../../user-guide/reference-data.md>); [Source](<https://github.com/bulik/ldsc/issues/33#issuecomment-105040547>).

<a id="issue-45"></a>

### #45: HapMap3 SNPs link

[Original issue](https://github.com/bulik/ldsc/issues/45) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** An old /~bulik/w_hm3.snplist.bz2 URL returns a permission error.
Users cannot perform the default preprocessing tutorial.

**Legacy disposition.** Maintainer supplied the alkesgroup/LDSCORE replacement and states the link was corrected; this is historical evidence, not a current availability check.

**LDSC3 assessment — Resource-only concern.** LDSC3 bundles its HM3 map, so its default restriction does not depend on this particular download.
External legacy tutorials and custom workflows may still reference the obsolete host.
Bundled data require provenance and licensing documentation; packaging a resource does not establish that any similarly named online list is equivalent.

**Recommended prevention.** Mark the original URL historical and point users to the bundled default plus currently checked authoritative sources in the resource audit.

**Evidence and verification.** source-inspected
[src/ldsc/data/hm3_curated_map.tsv.gz](<../../../../src/ldsc/data/hm3_curated_map.tsv.gz>); [src/ldsc/data/ATTRIBUTION.txt](<../../../../src/ldsc/data/ATTRIBUTION.txt>); [docs/user-guide/reference-data.md](<../../../user-guide/reference-data.md>).

<a id="issue-46"></a>

### #46: 1000G.mac5eur.* files

[Original issue](https://github.com/bulik/ldsc/issues/46) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **3**.

Workflow tags: `ldscore`, `baseline-partitioned-h2`.

**Trigger and harm.** Tutorial mentions 1000G.mac5eur files but only tarball names appear at the download location.
Users mistake packaging for missing data or choose an incompatible archive.

**Legacy disposition.** Maintainer identifies 1000G_plinkfiles.tgz and 1000G_frq.tgz; later users still report discoverability problems.

**LDSC3 assessment — Resource-only concern.** LDSC3 can use independently prepared PLINK panels and does not need these historical filenames.
Reference-data guidance should identify required roles, archive members and genome build instead of assuming readers infer contents from archive names.
The same release must supply a coherent genotype, frequency/count and annotation set where required.

**Recommended prevention.** For each recommended resource provide a short archive-to-member map and compatibility note, without downloading large panels during this audit.

**Evidence and verification.** source-inspected
[docs/user-guide/reference-data.md](<../../../user-guide/reference-data.md>); [docs/reference/formats/legacy-ldscores.md](<../../../reference/formats/legacy-ldscores.md>).

<a id="issue-61"></a>

### #61: missing file from pat. heritability tutorial

[Original issue](https://github.com/bulik/ldsc/issues/61) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **3**.

Workflow tags: `baseline-partitioned-h2`.

**Trigger and harm.** The tutorial names members unavailable as individual links.
Users cannot select or extract required reference files confidently.

**Legacy disposition.** Users locate the members in tarballs and maintainer explicitly updates the wiki.

**LDSC3 assessment — Resource-only concern.** The historical files are not mandatory LDSC3 filenames.
The resource wiki should still explain archive members and companion roles, particularly annotations, counts, regression weights and frequencies.
Link checking alone will not catch this discoverability problem because the tarball URL can succeed while instructions remain unusable.

**Recommended prevention.** List each required role next to its archive/member name and a short inspect command; deduplicate this recommendation with #46 while retaining both issue records.

**Evidence and verification.** source-inspected
[docs/user-guide/reference-data.md](<../../../user-guide/reference-data.md>); [docs/reference/formats/legacy-ldscores.md](<../../../reference/formats/legacy-ldscores.md>); [Source](<https://github.com/bulik/ldsc/issues/61#issuecomment-256978663>).

<a id="issue-106"></a>

### #106: Test enrichment in new annotations

[Original issue](https://github.com/bulik/ldsc/issues/106) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **3**.

Workflow tags: `ldscore`, `ldsc-seg`.

**Trigger and harm.** Users cannot find named genotype, HM3 and CNS annotation members at the directory URL.
Beginners mix Phase 1 and Phase 3 resources or infer wrong files from old prefixes.

**Legacy disposition.** Maintainer supplies archive names; user notes actual prefix/release inconsistencies and ultimately reports success.

**LDSC3 assessment — Resource-only concern.** LDSC3 accepts independent coherent resources and uses declared source roles.
The resource wiki should map roles to exact release-specific archive contents, not reproduce obsolete 1000G.mac5eur naming as a contract.
Matching only an archive label is insufficient when members changed prefixes.

**Recommended prevention.** Provide release-specific example member names and compatibility groups, with a maintained access check date.

**Evidence and verification.** source-inspected
[docs/user-guide/reference-data.md](<../../../user-guide/reference-data.md>); [docs/reference/input-paths.md](<../../../reference/input-paths.md>); [docs/reference/formats/legacy-ldscores.md](<../../../reference/formats/legacy-ldscores.md>).

<a id="issue-169"></a>

### #169: local version vs LDhub

[Original issue](https://github.com/bulik/ldsc/issues/169) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `h2`.

**Trigger and harm.** Reference and retained SNP counts differ across services and local runs.
Users attribute differences to software when resources, QC or versions differ.

**Legacy disposition.** Open without an authoritative service-version answer.

**LDSC3 assessment — Resource-only concern.** LDSC3 stores resource/configuration metadata and retained counts, but cannot recover an external service snapshot from a trait label.
Matching panel release, HM3 definition, region filters, sample-size policy and estimator settings is necessary before expecting equality. A newer panel is not automatically the same estimand or the correct replication resource.

**Recommended prevention.** Provide a compact provenance comparison template and mark unavailable external-service versions explicitly rather than speculate which one is latest.

**Evidence and verification.** source-inspected
[docs/user-guide/reference-data.md](<../../../user-guide/reference-data.md>); [docs/reference/logs-and-outputs.md](<../../../reference/logs-and-outputs.md>); [src/ldsc/_regression/runner.py:_rg_pair_metadata](<../../../../src/ldsc/_regression/runner.py>).

<a id="issue-179"></a>

### #179: Tutorial data no longer available

[Original issue](https://github.com/bulik/ldsc/issues/179) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `munge-sumstats`, `h2`, `rg`.

**Trigger and harm.** Old med.unc.edu direct links for cross-disorder bipolar and schizophrenia examples return404. A novice cannot reproduce the first tutorial even if package installation succeeds.

**Legacy disposition.** A commenter supplies the PGC results portal and advises accepting dataset terms; that is a resource workaround, not a code fix.

**LDSC3 assessment — Resource-only concern.** LDSC3 needs suitable public summary statistics for demonstrations, not necessarily these historical filenames.
Its guided workflow cannot make a withdrawn direct link valid, and substitute studies require new provenance and sample-size settings.

**Recommended prevention.** Use a currently verified official landing page, dataset identifiers and access instructions; distinguish gated from dead links.

**Evidence and verification.** source-inspected
[Source](<https://github.com/bulik/ldsc/issues/179#issuecomment-556287565>); [docs/user-guide/guided-tutorial.md](<../../../user-guide/guided-tutorial.md>).

<a id="issue-201"></a>

### #201: Partitioned h2: Difference between GTEx_brain and Multi_tissue_gene_expr datasets

[Original issue](https://github.com/bulik/ldsc/issues/201) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** Brain amygdala LD scores differ between two similarly named distributed annotation collections.
Users may compare scores whose specificity definitions and control universes differ.

**Legacy disposition.** Reporter resolves the question through the2018 paper Methods subsection Analysis of13brain regions using data from GTEx; no code fix.

**LDSC3 assessment — Resource-only concern.** LDSC3 can consume either collection after compatible conversion, but matching tissue names do not establish identical query annotations.
The statistical meaning comes from the original gene ranking, expression comparison and control construction, which file names alone do not retain.

**Recommended prevention.** Record annotation-generation methods, tissue comparison universe and control definition in the later resource catalog.

**Evidence and verification.** source-inspected
[Source](<https://github.com/bulik/ldsc/issues/201#issuecomment-608026518>); [docs/user-guide/workflows/cell-type-enrichment.md](<../../../user-guide/workflows/cell-type-enrichment.md>).

<a id="issue-208"></a>

### #208: rg tutorial broken

[Original issue](https://github.com/bulik/ldsc/issues/208) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `rg`.

**Trigger and harm.** PGC direct BIP/SCZ downloads no longer work.
Tutorial onboarding fails before the statistical workflow can be learned.

**Legacy disposition.** Commenter directs users to issue179 and the PGC portal; no code change identified.

**LDSC3 assessment — Resource-only concern.** LDSC3 requires compatible demonstration summary statistics, not those exact historical archives.
Source access conditions and current landing pages must be documented; changing a link does not itself verify dataset equivalence or the tutorial N values.

**Recommended prevention.** Replace stale direct downloads with verified official landing pages and a clearly identified tutorial dataset.

**Evidence and verification.** source-inspected
[Source](<https://github.com/bulik/ldsc/issues/208>); [docs/user-guide/guided-tutorial.md](<../../../user-guide/guided-tutorial.md>).

<a id="issue-224"></a>

### #224: Manifest file for partitioned heritability

[Original issue](https://github.com/bulik/ldsc/issues/224) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `baseline-partitioned-h2`, `continuous`.

**Trigger and harm.** User asks what names such as MAF_Adj_LLD_AFR mean.
Opaque annotations impede interpretation and lead to wrong biological conclusions.

**Legacy disposition.** Open; no explanatory manifest is supplied in the thread.

**LDSC3 assessment — Resource-only concern.** LDSC3 preserves annotation column names and type metadata but cannot reconstruct their scientific provenance from abbreviations.
In particular AFR in this baseline annotation refers to the original annotation definition and should not be casually interpreted as the ancestry of the whole European reference panel.

**Recommended prevention.** Catalog original annotation definitions, paper/table references, transformations and ancestry roles; do not invent expanded meanings from names.

**Evidence and verification.** source-inspected
[docs/reference/formats/column-schema.md](<../../../reference/formats/column-schema.md>); [docs/user-guide/workflows/continuous-annotations.md](<../../../user-guide/workflows/continuous-annotations.md>).

<a id="issue-231"></a>

### #231: minor wget issue with the heritability example on the wiki

[Original issue](https://github.com/bulik/ldsc/issues/231) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** Two PGC tutorial wget commands return404 while LD-reference downloads work.
Users blame their local configuration for a withdrawn external dataset link.

**Legacy disposition.** Comment links to issue179; no package fix.

**LDSC3 assessment — Resource-only concern.** The problem is the tutorial summary-statistics distribution route, not the LD-score files or LDSC3 regression logic. A maintained example needs a stable official landing page and explicit access conditions, with alternative datasets labeled as alternatives rather than identical archives.

**Recommended prevention.** Link verified PGC dataset landing pages and explain terms-controlled access; maintain a small freely accessible smoke example.

**Evidence and verification.** source-inspected
[Source](<https://github.com/bulik/ldsc/issues/231>); [docs/user-guide/guided-tutorial.md](<../../../user-guide/guided-tutorial.md>).

<a id="issue-273"></a>

### #273: Heritability tutorial download links broken

[Original issue](https://github.com/bulik/ldsc/issues/273) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** Historical bipolar/schizophrenia wget links return404 after the files move behind access terms.
Tutorials cannot be reproduced by copying the commands.

**Legacy disposition.** Comment links an official PGC form and suggests a different example dataset; no code fix.

**LDSC3 assessment — Resource-only concern.** LDSC3 should document access requirements and dataset identity rather than promise unattended direct download. A replacement pair should have openly accessible provenance and scientifically correct N/effect coding; it need not reproduce identical historical estimates.

**Recommended prevention.** Prefer stable official landing pages and a small freely accessible tutorial pair with fixed expected behavior.

**Evidence and verification.** source-inspected
[Source](<https://github.com/bulik/ldsc/issues/273#issuecomment-783875791>); [docs/user-guide/guided-tutorial.md](<../../../user-guide/guided-tutorial.md>).

<a id="issue-282"></a>

### #282: Data for tutorial no longer available

[Original issue](https://github.com/bulik/ldsc/issues/282) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** PGC direct download404 is followed by a suggested gist that later users cannot use.
An unverified replacement link prolongs tutorial failure.

**Legacy disposition.** No permanent resolution; the purported mirror is challenged in a later comment.

**LDSC3 assessment — Resource-only concern.** LDSC3 tutorial resources should be identified by official study/accession rather than relying on a historical gist.
Current availability, reuse conditions and correspondence to the tutorial must each be checked; a comment saying a mirror exists is only a lead.

**Recommended prevention.** Use verified official landing pages and maintain an openly accessible example with explicit provenance.

**Evidence and verification.** source-inspected
[Source](<https://github.com/bulik/ldsc/issues/282>); [docs/user-guide/guided-tutorial.md](<../../../user-guide/guided-tutorial.md>).

<a id="issue-311"></a>

### #311: LD Hub website issues

[Original issue](https://github.com/bulik/ldsc/issues/311) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** Users cannot access the LD Hub website or obtain a response from its operators. A hosted analysis workflow becomes inaccessible independently of local software functionality.

**Legacy disposition.** Open with corroborating users; no public service restoration or shutdown confirmation.

**LDSC3 assessment — Resource-only concern.** LDSC3 is a local package and can run analyses from available summary statistics and LD-score directories without LD Hub.
It does not recreate LD Hub hosted GWAS inventory, account/upload service or support operations, so website downtime should not be classified as a persistent local code defect.

**Recommended prevention.** Document the local alternative and required resources; record current service availability separately from historical outage reports.

**Evidence and verification.** source-inspected
[docs/user-guide/guided-tutorial.md](<../../../user-guide/guided-tutorial.md>); [docs/user-guide/workflows/genetic-correlation.md](<../../../user-guide/workflows/genetic-correlation.md>); [Source](<https://github.com/bulik/ldsc/issues/311>).

<a id="issue-325"></a>

### #325: LD Hub

[Original issue](https://github.com/bulik/ldsc/issues/325) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **4**.

Workflow tags: `rg`.

**Trigger and harm.** LD Hub remains inaccessible and users seek a multi-trait correlation alternative.
Dependence on a hosted service prevents analyses despite available local methods.

**Legacy disposition.** Comments suggest another website, but no restoration or scientific validation of an alternative is demonstrated.

**LDSC3 assessment — Resource-only concern.** LDSC3 offers local all-pairs genetic correlation and plotting, addressing the analysis mechanics when valid inputs are available.
It does not include the external hosted summary-statistics database, service guarantees or automatic equivalence to another portal.

**Recommended prevention.** Explain the local workflow and its resource requirements; verify any recommended portal and separate availability from analytical equivalence.

**Evidence and verification.** source-inspected
[docs/user-guide/workflows/genetic-correlation.md](<../../../user-guide/workflows/genetic-correlation.md>); [Source](<https://github.com/bulik/ldsc/issues/325>).

<a id="issue-330"></a>

### #330: HTTP request sent, awaiting response... 404 Not Found

[Original issue](https://github.com/bulik/ldsc/issues/330) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **4**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** HTTP redirects eventually lead to a removed historical bipolar/schizophrenia archive.
Users waste time changing network settings for an absent resource.

**Legacy disposition.** A Figshare mirror is suggested then reported invalid; no verified durable replacement.

**LDSC3 assessment — Resource-only concern.** LDSC3 should treat mirrors in comments as discovery leads and use a currently verified official study landing page. A replacement tutorial can use other datasets with explicit provenance; copying the old expected values or N would be misleading if the study changed.

**Recommended prevention.** Record final redirect/status and access conditions for each source; replace the tutorial dataset only with a fully identified compatible example.

**Evidence and verification.** source-inspected
[Source](<https://github.com/bulik/ldsc/issues/330>); [docs/user-guide/guided-tutorial.md](<../../../user-guide/guided-tutorial.md>).

<a id="issue-338"></a>

### #338: Heritability and Genetic Correlation wget Error 404

[Original issue](https://github.com/bulik/ldsc/issues/338) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** A redirected historical summary-statistics wget command fails with404.
The first tutorial cannot be followed from the documented command.

**Legacy disposition.** Closed as duplicate of issue330, not because the resource was restored.

**LDSC3 assessment — Resource-only concern.** LDSC3 should deduplicate this failure in recommendations while retaining the individual issue and its frequency contribution.
The closed duplicate state supplies no evidence of a working current URL or access permission; those require a live resource check.

**Recommended prevention.** Maintain one verified replacement/access explanation for the shared PGC failure and link every duplicate issue to it.

**Evidence and verification.** source-inspected
[Source](<https://github.com/bulik/ldsc/issues/338#issuecomment-1033143499>); [docs/user-guide/guided-tutorial.md](<../../../user-guide/guided-tutorial.md>).

<a id="issue-370"></a>

### #370: questions about pre-computed ldscore

[Original issue](https://github.com/bulik/ldsc/issues/370) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **6**.

Workflow tags: `ldscore`, `h2`.

**Trigger and harm.** Broad LDSCORE links return forbidden or fail during a hosting disruption.
Analyses cannot acquire their LD-score and HM3 inputs.

**Legacy disposition.** One user reports temporary recovery; others still fail and cite a cloud-cost discussion.
That is hosting evidence, not a code fix.

**LDSC3 assessment — Resource-only concern.** LDSC3 packages its HM3 map but still requires suitable LD resources. A reachable mirror needs provenance and compatibility checks, and 403 must not be described as definitive deletion.

**Recommended prevention.** Maintain a dated resource index with authoritative locations, availability labels and compatible bundles; retain the historical failure context.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/data/ATTRIBUTION.txt](<../../../../src/ldsc/data/ATTRIBUTION.txt>); [docs/reference/formats/legacy-ldscores.md](<../../../reference/formats/legacy-ldscores.md>); [Source](<https://github.com/bulik/ldsc/issues/370>).

<a id="issue-371"></a>

### #371: link to reference files does not work

[Original issue](https://github.com/bulik/ldsc/issues/371) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **11**.

Workflow tags: `munge-sumstats`, `h2`.

**Trigger and harm.** Users cannot find EUR scores and HM3 files after the Broad hosting move.
Missing migration instructions repeatedly block the tutorial and encourage unverified file sharing.

**Legacy disposition.** Comments point to a new-location notice and a Colorado teaching mirror; users confirm a download, not scientific equivalence.

**LDSC3 assessment — Resource-only concern.** The packaged HM3 map removes one external dependency, while legacy score conversion still needs complete paired suites and counts. A mirror should be identified as such rather than presented as an independently curated release.

**Recommended prevention.** List the role and provenance of each replacement resource and give direct bounded downloads rather than an unexplained recursive directory fetch.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_genome/hm3.py](<../../../../src/ldsc/_genome/hm3.py>); [src/ldsc/_ldscore/legacy_converter.py](<../../../../src/ldsc/_ldscore/legacy_converter.py>); [docs/reference/formats/legacy-ldscores.md](<../../../reference/formats/legacy-ldscores.md>).

<a id="issue-377"></a>

### #377: Reference files - 1000G LD scores

[Original issue](https://github.com/bulik/ldsc/issues/377) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `h2`, `munge-sumstats`.

**Trigger and harm.** Requests legacy EUR scores and HM3 list after the hosting URL fails.
The tutorial cannot be run without acquiring compatible inputs.

**Legacy disposition.** Open; later comments report the same blockage but supply no resolution.

**LDSC3 assessment — Resource-only concern.** LDSC3 includes an HM3 map while LD-score input remains external or user-computed.
The historical archive name is not itself a requirement; matching reference, regression SNPs, weights and counts are.

**Recommended prevention.** Provide resource roles and dated canonical download locations with alternatives clearly labeled as mirrors.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_genome/hm3.py](<../../../../src/ldsc/_genome/hm3.py>); [docs/reference/formats/legacy-ldscores.md](<../../../reference/formats/legacy-ldscores.md>).

<a id="issue-378"></a>

### #378: Unable to download from google bucket

[Original issue](https://github.com/bulik/ldsc/issues/378) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **6**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** Cahoy cell-type archive appears inaccessible in the cloud bucket.
Billing errors masquerade as absent scientific resources and trigger unnecessary mirror requests.

**Legacy disposition.** The commenter ultimately confirms the problem was their billing account and the bucket could be used; no package fix was needed.

**LDSC3 assessment — Resource-only concern.** LDSC3 does not manage Google billing and needs the resource only when reproducing that particular analysis. A suitable regenerated query model is possible, but must preserve its gene-control and baseline meaning.

**Recommended prevention.** Explain requester-pays access separately from stale links and list a free compatible alternative only when its provenance and contents are verified.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/user-guide/workflows/cell-type-enrichment.md](<../../../user-guide/workflows/cell-type-enrichment.md>); [Source](<https://github.com/bulik/ldsc/issues/378#issuecomment-1889826602>).

<a id="issue-382"></a>

### #382: eas_ldscores.tar.bz2

[Original issue](https://github.com/bulik/ldsc/issues/382) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **6**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** Users request the EAS LD-score archive after the original host is unavailable.
Informal email/Drive sharing can obscure release, integrity and matching weights.

**Legacy disposition.** Closed after file sharing; a comment links the requester-pays migration notice and another a community Drive copy.
No code fix is implicated.

**LDSC3 assessment — Resource-only concern.** LDSC3 still needs an ancestry-matched source or recomputed LD scores.
The fact that a shared archive worked for one user does not establish its build, panel identity or compatibility with current baselines.

**Recommended prevention.** Prefer a versioned archival record with build, ancestry, score/count/weight families; preserve community mirrors only as provenance-qualified fallbacks.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/reference/formats/legacy-ldscores.md](<../../../reference/formats/legacy-ldscores.md>); [src/ldsc/_ldscore/legacy_converter.py](<../../../../src/ldsc/_ldscore/legacy_converter.py>).

<a id="issue-402"></a>

### #402: Could anyone send me the "eas_ldscores.tar.bz2"? I wanted to calculate the ldsc of Asia, but I found I didn't have such file. And I can't download it from cloud. If anyone can send me it, my e-mail is chenqiwjc@gmail.com. Thank you so much.

[Original issue](https://github.com/bulik/ldsc/issues/402) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** Requests a missing EAS LD-score archive; a follow-up asks for a different EUR archive.
Similar filenames can hide different scientific roles and ancestry.

**Legacy disposition.** Open; a commenter suggests computing scores but supplies no validated bundle.

**LDSC3 assessment — Resource-only concern.** LDSC3 can compute scores from a suitable panel and import supported legacy suites, but an EAS score set is not interchangeable with a small EUR tutorial panel.

**Recommended prevention.** Separate EAS analysis resources from tutorial data, recording ancestry, build, chromosome scope, weights and counts for each.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_ldscore/calculator.py](<../../../../src/ldsc/_ldscore/calculator.py>); [src/ldsc/_ldscore/legacy_converter.py](<../../../../src/ldsc/_ldscore/legacy_converter.py>); [docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>).

<a id="issue-403"></a>

### #403: 1kg_eur.tar.bz2 not found

[Original issue](https://github.com/bulik/ldsc/issues/403) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **7**.

Workflow tags: `ldscore`.

**Trigger and harm.** The old 1kg_eur tutorial archive URL returns 404; a mirror later proves to contain only chr22.
Users may mistake a miniature demonstration panel for genome-wide scientific input.

**Legacy disposition.** Comments identify the hosting move and a GitHub mirror; another notes its chr22-only scope.
No code fix is involved.

**LDSC3 assessment — Resource-only concern.** LDSC3 supports small fixtures for demonstrations, while production index construction requires full autosomal coverage.
Resource scope must be explicit even when a download succeeds.

**Recommended prevention.** Label demonstration versus full-analysis resources in the wiki and show chromosome counts before running a genome-wide workflow.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[tests/fixtures/minimal_external_resources/README.md](<../../../../tests/fixtures/minimal_external_resources/README.md>); [docs/reference/gene-ldscore-index.md](<../../../reference/gene-ldscore-index.md>); [Source](<https://github.com/bulik/ldsc/issues/403#issuecomment-2828307753>).

<a id="issue-409"></a>

### #409: bucket seems to be unavailable

[Original issue](https://github.com/bulik/ldsc/issues/409) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `ldscore`.

**Trigger and harm.** Another report says the requester-pays bucket seems unavailable.
Users cannot distinguish object absence from access or billing setup.

**Legacy disposition.** Open; only a follow-up asking for a workaround.

**LDSC3 assessment — Resource-only concern.** The thread supplies no HTTP response or billing diagnostics, so it cannot establish that the resource is dead.
LDSC3 does not control cloud account access.

**Recommended prevention.** Record access-required/blocked separately from verified missing objects, and document payer-project setup only from current cloud documentation.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/reference/formats/legacy-ldscores.md](<../../../reference/formats/legacy-ldscores.md>); [Source](<https://github.com/bulik/ldsc/issues/409>).

<a id="issue-412"></a>

### #412: Cell type specific analyses

[Original issue](https://github.com/bulik/ldsc/issues/412) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** Cell-type score links at the original Broad location no longer work.
Users cannot reproduce the legacy Cahoy/SEG examples.

**Legacy disposition.** A commenter points to the requester-pays LDSC_SEG_ldscores directory; a later user still needs another route.

**LDSC3 assessment — Resource-only concern.** LDSC3 can build query annotations from gene lists/BED, but reproducing an old tutorial requires the same biological annotation definition and matching baseline/control. A replacement archive must satisfy those roles.

**Recommended prevention.** Provide a dated resource table with billing/access status and an alternative tiny local example that demonstrates the workflow without the entire historical archive.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/user-guide/workflows/cell-type-enrichment.md](<../../../user-guide/workflows/cell-type-enrichment.md>); [Source](<https://github.com/bulik/ldsc/issues/412#issuecomment-1899300845>).

<a id="issue-414"></a>

### #414: Need LDSC file

[Original issue](https://github.com/bulik/ldsc/issues/414) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `ldscore`.

**Trigger and harm.** Requests the full EUR Phase3 PLINK reference archive after a stale Broad link.
LD-score generation cannot start without a matched genotype panel.

**Legacy disposition.** A commenter points to Zenodo record 7768714; no software fix is needed.

**LDSC3 assessment — Resource-only concern.** LDSC3 accepts PLINK BED/BIM/FAM and verifies trios, but the record must still be checked for release, genome build, ancestry and matching frequencies/maps.

**Recommended prevention.** List genotype panel separately from precomputed scores, counts, maps and weights; verify file contents before describing a resource as a complete analysis bundle.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_runtime/paths.py:454](<../../../../src/ldsc/_runtime/paths.py#L454>); [docs/user-guide/workflows/ld-scores.md](<../../../user-guide/workflows/ld-scores.md>); [Source](<https://zenodo.org/records/7768714>).

<a id="issue-429"></a>

### #429: https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2 download problem?

[Original issue](https://github.com/bulik/ldsc/issues/429) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **4**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** A user obtains new LD-score files but does not know whether they match the old EUR archive, then sees different p-values.
Successful downloads can conceal changes in panel, weights, counts or implementation.

**Legacy disposition.** Comments link Zenodo 10515792 but admit uncertainty about equivalence; author confirms files downloaded, not matched scientific inputs.

**LDSC3 assessment — Resource-only concern.** LDSC3 conversion validates structure, but cannot establish historical equivalence from similar archive names.
The subsequent p-value change needs its own matched-input/version investigation.

**Recommended prevention.** Record exact resource releases and roles, then compare retained SNPs, estimates and SE under matched settings before attributing differences to software.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_ldscore/legacy_converter.py](<../../../../src/ldsc/_ldscore/legacy_converter.py>); [docs/reference/formats/legacy-ldscores.md](<../../../reference/formats/legacy-ldscores.md>); [docs/dev/validation/index.md](<../index.md>).

<a id="issue-450"></a>

### #450: Cannot run tutorial

[Original issue](https://github.com/bulik/ldsc/issues/450) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** The cell-type tutorial hardcodes a Broad archive URL that returns 404. A novice cannot complete the advertised example.

**Legacy disposition.** A commenter directs the user to the requester-pays bucket; no confirmation of successful access follows.

**LDSC3 assessment — Resource-only concern.** LDSC3 supports the workflow but still requires compatible biological inputs. A miniature self-contained example would reduce reliance on a paid historical archive without claiming it reproduces Cahoy results.

**Recommended prevention.** Provide dated canonical links, billing/access notes and a small local example; keep historical resource reproduction separate from basic installation validation.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/user-guide/workflows/cell-type-enrichment.md](<../../../user-guide/workflows/cell-type-enrichment.md>); [tests/fixtures/minimal_external_resources/README.md](<../../../../tests/fixtures/minimal_external_resources/README.md>).

<a id="issue-465"></a>

### #465: LD Hub Website Inaccessibility

[Original issue](https://github.com/bulik/ldsc/issues/465) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** The LD Hub web service cannot be reached.
Users lose a hosted analysis route and may conflate service availability with the LDSC method or local package.

**Legacy disposition.** Open; follow-ups confirm access problems without a verified service resolution.

**LDSC3 assessment — Resource-only concern.** LDSC3 runs locally and does not depend on LD Hub, but replacing a hosted trait catalogue requires its own data provenance and access strategy.

**Recommended prevention.** Separate hosted-service status from package installation and provide a local workflow using independently acquired permitted inputs; do not promise to restore the remote service.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[README.md](<../../../../README.md>); [docs/user-guide/guided-tutorial.md](<../../../user-guide/guided-tutorial.md>); [Source](<https://github.com/bulik/ldsc/issues/465>).

<a id="issue-499"></a>

### #499: East Asian LD scores file

[Original issue](https://github.com/bulik/ldsc/issues/499) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** Requests the exact eas_ldscores.tar.bz2 counterpart of an older EUR bundle and cannot find it in Zenodo.
The absence of one basename can be mistaken for absence of suitable EAS resources.

**Legacy disposition.** Open with no answer or verified mapping.

**LDSC3 assessment — Resource-only concern.** LDSC3 needs a compatible reference/weight/count tuple or a matched panel for recomputation, not a particular historical archive name. A multi-annotation EAS suite and separate weights may supply equivalent roles after validated conversion.

**Recommended prevention.** Publish an EAS role-to-file map with build, panel, score columns, weights and count conventions; distinguish verified combinations from merely listed archives.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/reference/formats/legacy-ldscores.md](<../../../reference/formats/legacy-ldscores.md>); [src/ldsc/_ldscore/legacy_converter.py](<../../../../src/ldsc/_ldscore/legacy_converter.py>); [Source](<https://zenodo.org/records/10515792>).

<a id="issue-219"></a>

### #219: Change of reference URL (Suggestion)

[Original issue](https://github.com/bulik/ldsc/issues/219) · Source state: **open** · Severity: **low** · Accessible comments reviewed: **0**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** The manuscript hyperlink leads to a collection rather than the LDSC paper.
Users cannot reliably locate the foundational methods source.

**Legacy disposition.** Reporter proposes the Nature article [Source](<https://www.nature.com/articles/ng.3211>); no code fix.

**LDSC3 assessment — Resource-only concern.** LDSC3 already keeps citation metadata, but each documentation reference should point to an exact DOI or publication landing page. A resource catalog should distinguish methodological references from analysis files and retain publication identity even if a hosting URL changes.

**Recommended prevention.** Use DOI or exact article landing-page links and verify their destinations during documentation review.

**Evidence and verification.** source-inspected
[CITATION.cff](<../../../../CITATION.cff>); [Source](<https://www.nature.com/articles/ng.3211>).

