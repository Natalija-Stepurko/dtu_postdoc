"""Stage 1a. Fetch the reference complex and report what is in it."""
import urllib.request, gzip, hashlib
from common import DATA, PDB_ID, save_json
import biotite.structure.io.pdbx as pdbx
import biotite.structure as struc

url = f"https://files.rcsb.org/download/{PDB_ID}.cif.gz"
raw = DATA / f"{PDB_ID}.cif"
if not raw.exists():
    with urllib.request.urlopen(url, timeout=60) as r:
        raw.write_bytes(gzip.decompress(r.read()))

sha = hashlib.sha256(raw.read_bytes()).hexdigest()[:16]
f = pdbx.CIFFile.read(raw)
block = f.block

arr = pdbx.get_structure(f, model=1, use_author_fields=True)
prot = arr[struc.filter_amino_acids(arr)]

meta = {
    "pdb_id": PDB_ID,
    "source_url": url,
    "sha256_16": sha,
    "title": block["struct"]["title"].as_item().strip(),
    "method": block["exptl"]["method"].as_item(),
    "resolution_A": float(block["refine"]["ls_d_res_high"].as_item()),
    "deposition_date": block["pdbx_database_status"]["recvd_initial_deposition_date"].as_item(),
}
chains = {}
for ch in sorted(set(prot.chain_id)):
    sub = prot[prot.chain_id == ch]
    res_ids = sorted(set(sub.res_id.tolist()))
    seq = struc.to_sequence(sub)[0][0] if hasattr(struc, "to_sequence") else ""
    chains[ch] = {"n_residues": len(res_ids),
                  "first_res": int(res_ids[0]), "last_res": int(res_ids[-1]),
                  "n_atoms": int(sub.array_length()),
                  "sequence": str(seq)}
meta["chains"] = chains
save_json("s01_structure.json", meta)
for k, v in meta.items():
    if k != "chains": print(f"  {k}: {v}")
for c, v in chains.items():
    print(f"  chain {c}: {v['n_residues']} res ({v['first_res']}-{v['last_res']}), {v['n_atoms']} atoms")
    print(f"    {v['sequence']}")
