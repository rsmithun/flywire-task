# FlyWire shared circuit search

## Result

The largest exact weakly connected directed induced subgraph I found across the provided connectomes is a 12-neuron circuit shared by MANC v1.2.1, MAOL v1.1, and MCNS v0.9.

The matched neurons are listed in [network.csv](network.csv). I validated that the induced directed edge sets are identical across all three datasets and that the selected subgraph is weakly connected.

## Technical approach

I treated the task as an exact common-induced-subgraph problem on unweighted directed graphs.

1. I streamed the five edge lists and computed node-overlap statistics.
2. The strongest candidate trio was MANC, MAOL, and MCNS, which share a large exact common-edge component.
3. I reduced the search to the largest weakly connected common-edge component shared by those three datasets.
4. I solved the exact connectivity-and-independence problem with a rooted mixed-integer program using SciPy/HiGHS.
5. I validated the final 12-node solution by re-reading the raw edge lists and checking that the induced directed subgraphs match exactly.

## Assumptions

- Edge weights were ignored; only the unweighted directed edge lists were used.
- The requested circuit was interpreted as a weakly connected directed induced subgraph.
- The dataset order in [network.csv](network.csv) is MANC, MAOL, MCNS.

## Reproduction

1. Download the five challenge edge lists into `data/` using these exact names:
	- `data/banc_626_edge_list.csv`
	- `data/fafb_783_edge_list.csv`
	- `data/manc_1.2.1_edge_list.csv`
	- `data/maol_1.1_edge_list.csv`
	- `data/mcns_0.9_edge_list.csv`
2. Install Python dependencies: `python -m pip install -r requirements.txt`
3. Recompute the exact circuit: `python scripts/reproduce_shared_circuit.py --data-root . --out network.csv`
4. Regenerate the figures: `python scripts/make_figures.py --data-root . --network network.csv --out figures`

The figure script downloads the public MANC annotations and skeletons automatically. Temporary downloads and caches stay under `.deps/` and `figures/.cache/` and are ignored by git.

## Repository contents

- [network.csv](network.csv) contains the matched neuron table.
- [science.md](science.md) contains the one-page scientific summary.
- [scripts/reproduce_shared_circuit.py](scripts/reproduce_shared_circuit.py) contains the analysis code.
- [scripts/make_figures.py](scripts/make_figures.py) contains the visualization code.
- [requirements.txt](requirements.txt) lists the Python dependencies.

## Data sources

- `data/banc_626_edge_list.csv`
- `data/fafb_783_edge_list.csv`
- `data/manc_1.2.1_edge_list.csv`
- `data/maol_1.1_edge_list.csv`
- `data/mcns_0.9_edge_list.csv`
- Public MANC annotations and skeletons used for the biological summary and figures