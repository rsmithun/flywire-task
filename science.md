# Shared tactile circuit across MANC, MAOL, and MCNS

## Circuit

I identified a 12-neuron directed induced subgraph that is exactly conserved across MANC v1.2.1, MAOL v1.1, and MCNS v0.9. The circuit is weakly connected and is centered on neuron `10114`, annotated in MANC as `vMS17`, a GABAergic bilateral interconnecting interneuron.

The strongest sensory entry point is `34613` (`SNta02`), a thoracic tactile bristle afferent. Downstream of that entry, the module contains a small set of ventral nerve cord intrinsic interneurons with mixed ipsilateral, contralateral, and bilateral morphology. The MANC annotations for the 12 shared neurons indicate a mixture of GABAergic, glutamatergic, and cholinergic cell types.

## Biological interpretation

The topology is consistent with a conserved tactile-processing and bilateral coordination microcircuit. The sensory afferent feeds a compact interneuron fan-out, while the GABAergic hub distributes to multiple local interneurons, including a reciprocal partner (`37836`). That reciprocal loop suggests local gating or stabilization rather than a purely feedforward path.

The MANC annotations strengthen this interpretation: `11595` and `21042` are marked as tactile interneurons, and several other nodes are classified as ventral-nerve-cord intrinsic interneurons with bilateral or ipsilateral restriction. Together, this points to a segmental integration module that may shape tactile responses before they are relayed into downstream motor or premotor pathways.

## Visuals

- Network graph: [figures/network.png](figures/network.png)
- MANC skeleton render: [figures/manc_skeletons_3d.png](figures/manc_skeletons_3d.png)

## Methods note

I ignored edge weights and used only the unweighted directed edge lists. The shared circuit was validated by comparing the induced directed edge sets over the selected 12 nodes and confirming exact equality across MANC, MAOL, and MCNS.

## Key annotations in MANC

| Root ID | Cell type | NT | Note |
| --- | --- | --- | --- |
| 34613 | SNta02 | acetylcholine | Thoracic bristle sensory afferent |
| 10114 | vMS17 | gaba | Bilateral interconnecting hub |
| 28795 | IN19A120 | glutamate | Bilateral restricted interneuron |
| 37836 | IN20A.22A092 | acetylcholine | Reciprocal partner |
| 21042 | IN23B065 | acetylcholine | Tactile interneuron |
| 11595 | IN06B016 | gaba | Tactile interneuron |

## Citations

- Takemura et al. (2024). *A Connectome of the Male Drosophila Ventral Nerve Cord*. eLife.
- Marin et al. (2024). *Systematic annotation of a complete adult male Drosophila nerve cord connectome*. eLife.
- Cheong et al. (2024). *Transforming descending input into behavior*. eLife.
- Eckstein et al. (2024). *Neurotransmitter classification from electron microscopy images at synaptic sites in Drosophila melanogaster*. Cell.
- Schlegel et al. (2024). *Whole-brain annotation and multi-connectome cell typing of Drosophila*. Nature.
- Berg et al. (2025). *Sexual dimorphism in the complete connectome of the Drosophila male central nervous system*. bioRxiv.