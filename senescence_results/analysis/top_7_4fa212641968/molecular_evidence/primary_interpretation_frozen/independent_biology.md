# Independent biological assessment of the fixed top seven

Primary assessment frozen on 2026-09-25, before opening any SenMayo result, coverage, membership or leading-edge file. The Reactome specificity appendix was examined separately after forming the raw-weight interpretations below. K=7 was fixed upstream from the activity-score curve, not from biological annotations. This note does not change the workflow.

Inputs: `analysis/top_7_4fa212641968/selected_programs.csv`, candidate activity/distribution and mouse-influence outputs; `analysis/screen/class_results.csv`; `molecular_evidence/gene_weights_specificity_and_ranks.csv.gz`, raw Reactome enrichment and leading-edge outputs, and the corresponding specificity Reactome outputs. Reactome q values below use the complete 952-by-7 family. Ranks refer to the common 22,671-gene enrichment universe. A plus sign means positive dictionary loading, not a separately established age effect for that gene. These are associations of a learned component, not evidence that its pathways were experimentally activated.

## Main independent conclusion

Program 188 is the most biologically plausible senescence-compatible candidate among the seven, despite being fifth by the age-tail score. Its combination of lysosomal/myeloid, lipid-handling and inflammatory genes with a positive cell-cycle inhibitor is more informative than any individual marker. Activated/lipid-processing microglia remain a strong competing explanation; neither irreversible arrest nor senescence is established. Program 122 is a weaker, unresolved candidate whose more direct interpretation is inflammatory microglial activation. The three oligodendrocyte and two vascular programs have convincing lineage/age-associated interpretations, but the available evidence does not specifically support senescence. This is a useful distinction produced by interpretation after marker-independent prioritization, not a failure requiring a different K.

There is no basis here for claiming a newly discovered cell state. The most defensible novelty is the unsupervised program representation and marker-independent prioritization of interpretable aging components, including a candidate resembling previously described inflammatory microglia. Correspondence with the source atlas is contextual consistency, not independent replication.

## Compact assessment and exact supporting evidence

Each listed support set contains 3–5 positive leading-edge genes from the two stated primary Reactome themes. Other strong lineage genes are reported separately as context, rather than silently presented as Reactome leading-edge support.

| Program (score rank) | Dominant class; normalized score | Two representative primary functional themes, with joint q | Strong leading-edge support: signed raw loading (rank) | Strongest alternative; frozen senescence judgment |
|---|---|---|---|---|
| 179 (1) | OPC–Oligo; 0.715910 | Prostaglandin/thromboxane synthesis, R-MMU-2162123, q=0.014987; axon guidance, R-MMU-422475, q=0.000425 | Ptgds +0.080387 (9); Tubb4a +0.052757 (54); Dnm3 +0.047246 (83); Nfasc +0.045131 (94) | Mature oligodendrocyte/axon-associated component with prostanoid metabolism; **senescence not specifically supported**. |
| 163 (2) | OPC–Oligo; 0.588256 | G-alpha-i signaling, R-MMU-418594, q=7.76e-7; iron uptake/transport, R-MMU-917937, q=0.042567 | Fth1 +0.094515 (1); Grm7 +0.085530 (4); Trf +0.068856 (10); App +0.058693 (25); Pde4b +0.056366 (30) | Iron-handling/receptor-signaling oligodendrocyte component or finer maturation identity; **senescence not specifically supported**. |
| 122 (3) | Immune; 0.556530 | Cytokine signaling, R-MMU-1280215, q=0.002174; Toll-like receptor cascades, R-MMU-168898, q=0.013914 | Socs3 +0.069426 (6); Fos +0.059206 (22); Ly86 +0.058925 (24); Jun +0.058485 (25) | Inflammatory microglial activation with feedback regulation; **unresolved, weaker senescence plausibility than 188**. |
| 196 (4) | OPC–Oligo; 0.473634 | RHO GTPase cycle, R-MMU-9012999, q=1.83e-7; axon guidance, R-MMU-422475, q=0.000677 | Dpysl2 +0.048249 (39); Tubb4a +0.045248 (48); Dock10 +0.042110 (62); Cdc42bpa +0.041736 (66) | Mature oligodendrocyte structural/cytoskeletal component; **senescence not specifically supported**. |
| 188 (5) | Immune; 0.412978 | Lysosomal/myeloid machinery represented by neutrophil degranulation, R-MMU-6798695, q=5.70e-9; TP53 transcriptional regulation, R-MMU-3700989, q=0.000425 | Ctsd +0.105065 (1); Ctsb +0.096244 (2); Tyrobp +0.056445 (25); Cdkn1a +0.037350 (153) | Activated, lipid-processing inflammatory microglia; **plausible senescence-compatible candidate, not confirmed senescence**. |
| 60 (6) | Vascular; 0.365839 | VEGF signaling, R-MMU-194138, q=0.000161; WNT signaling, R-MMU-195721, q=0.000967 | Flt1 +0.083277 (1); Lef1 +0.062413 (12); Plcb1 +0.055204 (24); Kdr +0.049293 (39) | Endothelial growth-factor/barrier signaling associated with aging; **senescence not specifically supported**. |
| 124 (7) | Vascular; 0.316102 | Muscle contraction, R-MMU-397014, q=5.70e-9; potassium channels, R-MMU-1296071, q=0.005529 | Abcc9 +0.092259 (4); Cald1 +0.081099 (8); Kcnj8 +0.065411 (19); Myl9 +0.059154 (30) | Pericyte/mural-cell contractile and ion-channel identity; **senescence not specifically supported**. |

The term “neutrophil degranulation” should be explained as shared lysosomal/myeloid machinery in program 188. It is not evidence that this immune-class component represents neutrophils. Likewise, TP53 enrichment is not proof of irreversible cell-cycle arrest.

### Program-specific context that changes interpretation

**179, 163 and 196:** All contain strongly positive myelin/oligodendrocyte genes: Plp1 ranks 2 in every component; Mbp ranks 10, 15 and 5, respectively. Their functional subdivisions should be described as components within the lineage, not as three proven distinct senescent populations. For 179, Ptgds is prominent, but prostaglandin-pathway support contains only two leading genes (Ptgds and Ptges3), so avoid an extensive prostanoid mechanism story. For 163, Fth1/Trf give a concrete iron-handling interpretation, but the pathway q=0.0426 is modest and ranking-sensitive. For 196, Scd2, Mal, Mobp, Apod and Il33 provide additional mature-glial context. Il33 alone does not establish a pathological secretory phenotype: constitutive oligodendrocyte IL-33 has a demonstrated physiological/injury-response role [Gadani et al., 2015](https://doi.org/10.1016/j.neuron.2015.01.013).

Program 163 has a negative Opalin loading (−0.018230, rank 22,670). It would be tempting to equate this with the source atlas's Opalin-low old oligodendrocyte cluster, but that is not warranted: Cdkn1a and Art3 are negative here. Across 179/163/196, Cdkn1a is weak/negative rather than a prominent common support gene. Negative loading is a feature of the fitted component, not a direct age-related decrease. Broad-class adjustment cannot distinguish a shifted mature-oligodendrocyte mixture from altered activity within a finer subtype; one general limitation sentence is sufficient.

**122 versus 188:** The two immune programs are not interchangeable. Program 122 combines Socs3/AP-1-associated genes with strong P2ry12 (+0.065498, rank 10), Cx3cr1 (+0.063223, rank 13), and Siglech (+0.065337, rank 11). Ccl4 is also strong (+0.066818, rank 9), while Cdkn1a is less prominent (+0.015100, rank 835). Thus inflammatory signaling with retained microglial identity is the cleaner explanation; no need to invent senescence from generic stress enrichment. A Ccl4-positive inflammatory microglial subset can expand with aging or injury [Hammond et al., 2019](https://pubmed.ncbi.nlm.nih.gov/30471926/).

Program 188 has a substantially more coherent senescence-compatible combination. Beyond the four strictly compliant table genes, Apoe (+0.081794, rank 6), Lpl (+0.080581, rank 7), Tnf (+0.064395, rank 14), Ccl3 (+0.061418, rank 18), Bcl2a1d (+0.037992, rank 143), Spp1 (+0.037675, rank 149), Cst7 (+0.033683, rank 198), Itgax (+0.033062, rank 204) and Il1b (+0.026441, rank 353) contextualize lipid processing, inflammation and survival. Supporting secondary Reactome results are cytokine signaling (q=0.000139), Toll-like receptor cascades (q=0.001260), and lipoprotein assembly/remodeling/clearance (q=0.005559). These can substantiate a short mechanistic paragraph without adding more columns or treating every overlapping pathway as another independent vote. The inflammatory activation alternative remains credible even with Cdkn1a and survival genes; the thesis's existing Hall et al. citation is useful for the broader point that senescence-associated myeloid features can be reversible activation phenomena.

**60 versus 124:** Program 60 is endothelial-like, with Flt1, Igf1r (+0.079415, rank 2), Hdac9 (+0.072777, rank 3), Cxcl12 (+0.069526, rank 5), Ptprb (+0.062532, rank 11) and Lef1. A primary aging-endothelium study independently found age-associated Flt1/Cxcl12 and growth-factor/immune changes, and demonstrated responsiveness to circulatory cues [Chen et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7292569/). The current component supports that endothelial context, not demonstrated senescent arrest. Program 124's Rgs5 (+0.100948, rank 2), Abcc9, Pdgfrb (+0.076442, rank 10), Kcnj8 and Notch3 (+0.064014, rank 21) distinguish mural/pericyte-like identity from endothelial program 60. The vascular atlas supports treating these as different vascular populations [Vanlandewijck et al., 2018](https://www.nature.com/articles/nature25739). Because the regression adjusts only the broad vascular class, an endothelial/pericyte mixture shift is still a possible contributor.

## Relationship to the original atlas, without a novelty overclaim

The source atlas already reports age-enriched inflammatory microglia with Lpl/Cst7 and Cdkn1a/Bcl2a1 genes, increased Hdac9 in aged endothelial cells, and Dpyd/Abca8a changes in mature oligodendrocytes. Its inflammatory microglial cluster 843 is related to Hammond's OA2 population. Program 188 resembles that known pattern; 60 and the oligodendrocyte programs also overlap known aging biology. Conversely, 163 should not be identified with the atlas's Opalin-low, Art3/Cdkn1a-positive oligodendrocyte cluster solely from its negative Opalin loading. These are gene-level contextual similarities, not a tested mapping between learned programs and atlas clusters. [Jin et al., 2025](https://www.nature.com/articles/s41586-024-08350-8).

## Activity evidence and its limits

All seven are class-restricted under the agreed identity/activity rule. Dominant-class usage weights are 0.760, 0.894, 0.805, 0.842, 0.637, 0.909 and 0.960 in score order. Nearly the whole positive score comes from each dominant class: its normalized contributions are 0.709444, 0.588294, 0.555534, 0.468611, 0.402161, 0.364937 and 0.316119. Immune support therefore should not be described as an atlas-wide senescence program. These dominant classes have substantial donor coverage (OPC–Oligo 53 young/37 aged; Immune 26/27; Vascular 27/34), rather than the three-mice-per-age classes that were a concern in planning.

For each candidate, the fitted upper-tail age increase exceeds the positive overall-mean age increase. The dominant-class (mean, tail) coefficients are respectively: 179 (0.020401, 0.081705), 163 (0.017472, 0.070281), 122 (0.016366, 0.056798), 196 (0.016412, 0.050764), 188 (0.011292, 0.042852), 60 (0.033042, 0.084107), 124 (0.035597, 0.112295). They support the planned descriptive upper-tail criterion, not specificity for senescence.

Every single-mouse omission preserves positive scores. Ranks 1–6 are unchanged; 124 varies between ranks 7 and 8. Score ranges: 179 [0.703348, 0.726312], 163 [0.572008, 0.601842], 122 [0.528824, 0.577666], 196 [0.455727, 0.489735], 188 [0.370052, 0.441205], 60 [0.356945, 0.377788], 124 [0.305353, 0.329296]. These are influence ranges, not confidence intervals or sampling uncertainty.

The standardized distribution's zero fractions fall from young to aged for every candidate. For 122 they fall from 0.5876 to 0.4400; for 60 from 0.6511 to 0.4563. Thus many cells can have nonzero activity even though the score measures an upper-15% excess. Do not describe nonzero cells as senescent cells or claim that exactly 15% of cells form a discovered senescent population.

## Reactome specificity appendix, assessed before SenMayo

The appendix asks whether the same primary themes remain supported under gene-column-normalized weights; it is not a second discovery screen or independent replication. Raw and specificity q values belong to separate prespecified complete families. Do not choose whichever is smaller.

| Program | What survives for the prechosen primary interpretation | What weakens |
|---|---|---|
| 179 | Axon guidance q=0.0281; RHO cycle q=0.000587 | Prostaglandin synthesis q=0.1385, so that detail is ranking-sensitive. |
| 163 | Broad lineage interpretation remains in raw genes; potassium-channel q=0.0378 and RHO q=0.0445 remain as context | Primary G-alpha-i q=0.0949 and iron transport q=0.346 lose adjusted support. Do not call its specific iron/GPCR interpretation robust. |
| 122 | Cytokine q=0.00130 and TLR q=0.00527 | No reversal of the primary activation interpretation; no basis to upgrade it to senescence. |
| 196 | RHO cycle q=0.000450 | Axon guidance q=0.1415; keep the narrower cytoskeletal/lineage description. |
| 188 | Myeloid/lysosomal q=1.42e-8; TP53 q=0.00291. Supporting cytokine q=2.19e-6, TLR q=2.32e-6 and lipoprotein q=0.00504 also remain | The specific iron-transport detail weakens (q=0.0690); the overall mechanistic interpretation survives. |
| 60 | VEGF q=0.00159 | WNT q=0.1385; endothelial signaling is supported more strongly than a specific WNT mechanism. |
| 124 | Contraction q=1.42e-8 and potassium channels q=2.15e-6 | No material weakening of its primary contractile/mural interpretation. |

Translation, ribosomal and RNA-processing pathways dominate many primary and appendix results, even after column normalization. These annotations are real outputs but are not a sufficient discriminating senescence mechanism. None of the seven has q<=0.05 in the directly labeled Reactome Cellular Senescence branch in the raw analysis; the parent pathway is also unsupported in the specificity analysis. This does not veto candidate 188, but it prevents describing Reactome as a direct senescence validation.

## Citation budget proposal

Jin et al. and the general senescence/activation references are already in the thesis. Add at most these four biology papers, and omit Gadani if Il33 does not enter the main narrative:

1. Hammond et al. (2019), *Immunity* 50:253–271.e6, DOI 10.1016/j.immuni.2018.11.004 — inflammatory microglial aging/injury alternative for 122 and 188.
2. Chen et al. (2020), *Cell Reports* 30:4418–4432.e4, DOI 10.1016/j.celrep.2020.03.012 — endothelial aging context for 60. First author is Michelle B. Chen, not Yang.
3. Vanlandewijck et al. (2018), *Nature* 554:475–480, DOI 10.1038/nature25739 — endothelial versus mural/pericyte interpretation for 60/124.
4. Gadani et al. (2015), *Neuron* 85:703–709, DOI 10.1016/j.neuron.2015.01.013 — constitutive oligodendrocyte IL-33 is not intrinsically a senescence signature.

No workflow-invalidating finding was identified. The main writing risk is overinterpreting generic aging/lineage programs or treating known inflammatory microglia as a novel, confirmed senescent cell state.
